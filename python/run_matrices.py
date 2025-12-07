#!/usr/bin/env python3
"""
Script to run all matrices from matrices.json using the FRACESSA Python bindings
and save results to a JSON file. Processes matrices sequentially.
"""

import json
import sys
import time
import csv
from datetime import datetime
from pathlib import Path
from fracessa_py import Fracessa, Matrix, FracessaError


def candidate_to_comparison_key(candidate):
    """
    Convert a candidate dict to a hashable tuple for comparison.
    Excludes candidate_id since it may differ between runs.
    Excludes reason_ess/stability since it may differ between implementations.
    """
    # Normalize vector to string for comparison
    vector = candidate.get('vector', '')
    if isinstance(vector, list):
        vector = ','.join(str(v) for v in vector)
    
    return (
        str(vector),
        str(candidate.get('support', '')),
        str(candidate.get('support_size', '')),
        str(candidate.get('extended_support', '')),
        str(candidate.get('extended_support_size', '')),
        str(candidate.get('shift_reference', '')),
        str(candidate.get('is_ess', '')),
        # reason_ess/stability excluded from comparison
        str(candidate.get('payoff', '')),
        str(candidate.get('payoff_double', ''))
    )


def compare_runtimes(results, baseline_file):
    """
    Compare runtimes of current run against baseline JSON file.
    
    Args:
        results: List of result dicts from current run
        baseline_file: Path to baseline JSON file
        
    Returns:
        dict with runtime comparison results
    """
    # Check if baseline file exists
    if not baseline_file.exists():
        return {
            "status": "SKIPPED",
            "reason": f"Baseline file not found: {baseline_file}",
            "compared": 0,
            "faster": 0,
            "slower": 0,
            "same": 0,
            "details": []
        }
    
    # Load baseline results
    try:
        with open(baseline_file, 'r') as f:
            baseline_data = json.load(f)
        baseline_results = baseline_data.get('results', [])
    except Exception as e:
        return {
            "status": "ERROR",
            "reason": f"Failed to load baseline file: {e}",
            "compared": 0,
            "faster": 0,
            "slower": 0,
            "same": 0,
            "details": []
        }
    
    # Create lookup dict for baseline timings
    baseline_timings = {}
    for result in baseline_results:
        matrix_id = result.get('id')
        if result.get('result', {}).get('success'):
            baseline_timings[matrix_id] = result['result'].get('timing', 0.0)
    
    # Compare current runtimes with baseline
    compared = 0
    faster = 0
    slower = 0
    same = 0
    details = []
    
    for result in results:
        matrix_id = result.get('id')
        current_result = result.get('result', {})
        
        if not current_result.get('success'):
            continue  # Skip failed matrices
        
        if matrix_id not in baseline_timings:
            continue  # Skip matrices not in baseline
        
        current_timing = current_result.get('timing', 0.0)
        baseline_timing = baseline_timings[matrix_id]
        
        compared += 1
        
        if baseline_timing == 0:
            # Avoid division by zero
            if current_timing == 0:
                same += 1
                speedup = 1.0
                percentage_change = 0.0
            else:
                slower += 1
                speedup = float('inf')
                percentage_change = float('inf')
        else:
            speedup = baseline_timing / current_timing if current_timing > 0 else float('inf')
            
            # Calculate percentage change: ((current - baseline) / baseline) * 100
            percentage_change = ((current_timing - baseline_timing) / baseline_timing) * 100
            
            if abs(current_timing - baseline_timing) < 0.001:  # Consider same if within 1ms
                same += 1
            elif current_timing < baseline_timing:
                faster += 1
            else:
                slower += 1
        
        details.append({
            "matrix_id": matrix_id,
            "dimension": result.get('dimension'),
            "current_timing": current_timing,
            "baseline_timing": baseline_timing,
            "speedup": speedup if speedup != float('inf') else None,
            "slower_by": current_timing - baseline_timing if current_timing > baseline_timing else None,
            "faster_by": baseline_timing - current_timing if current_timing < baseline_timing else None,
            "percentage_change": percentage_change if percentage_change != float('inf') else None
        })
    
    # Sort details by speedup ratio (best improvements first)
    details.sort(key=lambda x: x.get('speedup') or 0, reverse=True)
    
    return {
        "status": "SUCCESS",
        "compared": compared,
        "faster": faster,
        "slower": slower,
        "same": same,
        "details": details
    }


def compare_ess_counts(results, baseline_file):
    """
    Compare ESS counts of current run against baseline JSON file.
    
    Args:
        results: List of result dicts from current run
        baseline_file: Path to baseline JSON file
        
    Returns:
        dict with ESS count comparison results
    """
    # Check if baseline file exists
    if not baseline_file.exists():
        return {
            "status": "SKIPPED",
            "reason": f"Baseline file not found: {baseline_file}",
            "compared": 0,
            "mismatches": []
        }
    
    # Load baseline results
    try:
        with open(baseline_file, 'r') as f:
            baseline_data = json.load(f)
        baseline_results = baseline_data.get('results', [])
    except Exception as e:
        return {
            "status": "ERROR",
            "reason": f"Failed to load baseline file: {e}",
            "compared": 0,
            "mismatches": []
        }
    
    # Create lookup dict for baseline ESS counts
    baseline_ess_counts = {}
    for result in baseline_results:
        matrix_id = result.get('id')
        if result.get('result', {}).get('success'):
            baseline_ess_counts[matrix_id] = result['result'].get('ess_count', 0)
    
    # Compare current ESS counts with baseline
    compared = 0
    mismatches = []
    
    for result in results:
        matrix_id = result.get('id')
        current_result = result.get('result', {})
        
        if not current_result.get('success'):
            continue  # Skip failed matrices
        
        if matrix_id not in baseline_ess_counts:
            continue  # Skip matrices not in baseline
        
        current_ess_count = current_result.get('ess_count', 0)
        baseline_ess_count = baseline_ess_counts[matrix_id]
        
        compared += 1
        
        if current_ess_count != baseline_ess_count:
            mismatches.append({
                "matrix_id": matrix_id,
                "dimension": result.get('dimension'),
                "current_ess": current_ess_count,
                "baseline_ess": baseline_ess_count
            })
    
    return {
        "status": "SUCCESS",
        "compared": compared,
        "mismatches": mismatches
    }


def verify_against_baseline(all_candidates, baseline_file):
    """
    Verify produced candidates against baseline CSV.
    
    Args:
        all_candidates: List of candidate dicts from current run
        baseline_file: Path to baseline CSV file
        
    Returns:
        dict with verification results
    """
    from collections import defaultdict
    
    # Check if baseline file exists
    if not baseline_file.exists():
        return {
            "status": "SKIPPED",
            "reason": f"Baseline file not found: {baseline_file}",
            "verified": 0,
            "passed": 0,
            "failed": 0,
            "details": []
        }
    
    # Load baseline candidates
    baseline_candidates = []
    with open(baseline_file, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            baseline_candidates.append(row)
    
    # Group candidates by matrix_id
    produced_by_matrix = defaultdict(set)
    baseline_by_matrix = defaultdict(set)
    produced_candidate_counts = defaultdict(int)
    baseline_candidate_counts = defaultdict(int)
    produced_ess_counts = defaultdict(int)
    baseline_ess_counts = defaultdict(int)
    
    for candidate in all_candidates:
        matrix_id = candidate.get('matrix_id')
        if matrix_id is None:
            continue  # Skip candidates without matrix_id
        matrix_id = int(matrix_id)  # Ensure int type for comparison
        key = candidate_to_comparison_key(candidate)
        produced_by_matrix[matrix_id].add(key)
        produced_candidate_counts[matrix_id] += 1
        # Count ESS candidates
        if candidate.get('is_ess', ''):
            is_ess_str = str(candidate.get('is_ess', '')).lower()
            if is_ess_str in ('true', '1', 'yes'):
                produced_ess_counts[matrix_id] += 1
    
    for candidate in baseline_candidates:
        matrix_id = int(candidate.get('matrix_id'))
        key = candidate_to_comparison_key(candidate)
        baseline_by_matrix[matrix_id].add(key)
        baseline_candidate_counts[matrix_id] += 1
        # Count ESS candidates
        if candidate.get('is_ess', ''):
            is_ess_str = str(candidate.get('is_ess', '')).lower()
            if is_ess_str in ('true', '1', 'yes'):
                baseline_ess_counts[matrix_id] += 1
    
    # Categorize matrix_ids
    produced_ids = set(produced_by_matrix.keys())
    baseline_ids = set(baseline_by_matrix.keys())
    
    # Matrices in produced but not in baseline - ERROR
    not_in_baseline = produced_ids - baseline_ids
    # Matrices in baseline but not in produced - just skipped (not an error)
    skipped = baseline_ids - produced_ids
    # Matrices in both - need to verify
    to_verify = produced_ids & baseline_ids
    
    passed = 0
    failed = 0
    mismatches = []
    errors = []
    candidate_count_mismatches = []
    ess_count_mismatches = []
    
    # Check for matrices not in baseline (ERROR)
    for matrix_id in sorted(not_in_baseline):
        failed += 1
        errors.append({
            "matrix_id": matrix_id,
            "error": "Matrix ID not found in baseline"
        })
    
    # Verify matrices that are in both
    for matrix_id in sorted(to_verify):
        produced_set = produced_by_matrix[matrix_id]
        baseline_set = baseline_by_matrix[matrix_id]
        
        # Check candidate count
        produced_candidate_count = produced_candidate_counts[matrix_id]
        baseline_candidate_count = baseline_candidate_counts[matrix_id]
        if produced_candidate_count != baseline_candidate_count:
            failed += 1
            candidate_count_mismatches.append({
                "matrix_id": matrix_id,
                "produced_count": produced_candidate_count,
                "baseline_count": baseline_candidate_count
            })
        
        # Check ESS count
        produced_ess_count = produced_ess_counts[matrix_id]
        baseline_ess_count = baseline_ess_counts[matrix_id]
        if produced_ess_count != baseline_ess_count:
            failed += 1
            ess_count_mismatches.append({
                "matrix_id": matrix_id,
                "produced_ess_count": produced_ess_count,
                "baseline_ess_count": baseline_ess_count
            })
        
        # Check candidate sets (excluding reason_ess/stability)
        if produced_set == baseline_set:
            # Only count as passed if candidate count also matches
            if produced_candidate_count == baseline_candidate_count:
                passed += 1
        else:
            failed += 1
            extra = produced_set - baseline_set
            missing = baseline_set - produced_set
            mismatches.append({
                "matrix_id": matrix_id,
                "produced_count": len(produced_set),
                "baseline_count": len(baseline_set),
                "extra_count": len(extra),
                "missing_count": len(missing)
            })
    
    return {
        "status": "PASS" if failed == 0 else "FAIL",
        "verified": len(to_verify),
        "passed": passed,
        "failed": failed,
        "skipped": list(sorted(skipped)),
        "errors": errors,
        "mismatches": mismatches,
        "candidate_count_mismatches": candidate_count_mismatches,
        "ess_count_mismatches": ess_count_mismatches
    }


def process_matrix(matrix_data, executable_path):
    """
    Process a single matrix.
    
    Args:
        matrix_data: dict with matrix information
        executable_path: path to the fracessa executable
        
    Returns:
        tuple: (matrix_id, result_entry, candidates_list)
    """
    matrix_id = matrix_data['id']
    dimension = matrix_data['dimension']
    number_ess = matrix_data['number_ess']
    is_cs = matrix_data['is_cs']
    matrix_str = matrix_data['matrix']
    
    try:
        # Initialize FRACESSA interface
        fracessa = Fracessa(str(executable_path))
        
        # Create matrix object
        matrix = Matrix(matrix_str, dimension, is_circular=is_cs)
        
        # Run ESS computation
        result = fracessa.compute_ess(
            matrix=matrix,
            include_candidates=True,
            enable_logging=False,
            timeout=1800.0,  # 30 minutes timeout
            matrix_id=matrix_id
        )
        
        if result.success:
            # Convert candidates to the expected format
            candidates_data = []
            for candidate in result.candidates:
                candidate_dict = {
                    'candidate_id': candidate.candidate_id,
                    'vector': candidate.vector,
                    'support': candidate.support,
                    'support_size': candidate.support_size,
                    'extended_support': candidate.extended_support,
                    'extended_support_size': candidate.extended_support_size,
                    'shift_reference': candidate.shift_reference,
                    'is_ess': candidate.is_ess,
                    'reason_ess': str(candidate.reason_ess),
                    'payoff': candidate.payoff,
                    'payoff_double': candidate.payoff_double
                }
                candidates_data.append(candidate_dict)
            
            computation_result = {
                "success": True,
                "ess_count": result.ess_count,
                "timing": result.computation_time,
                "candidates": candidates_data
            }
        else:
            computation_result = {
                "success": False,
                "error": result.error,
                "timing": result.computation_time
            }
    
    except FracessaError as e:
        computation_result = {
            "success": False,
            "error": f"FRACESSA Error: {str(e)}",
            "timing": 0.0
        }
    except Exception as e:
        computation_result = {
            "success": False,
            "error": f"Exception: {str(e)}",
            "timing": 0.0
        }
    
    # Extract candidates for CSV (with matrix_id)
    candidates_list = []
    if computation_result["success"] and "candidates" in computation_result:
        for candidate in computation_result["candidates"]:
            candidate_with_matrix_id = {"matrix_id": matrix_id}
            candidate_with_matrix_id.update(candidate)
            candidates_list.append(candidate_with_matrix_id)
    
    # Build result entry without candidates for JSON
    result_without_candidates = {k: v for k, v in computation_result.items() if k != "candidates"}
    
    result_entry = {
        "id": matrix_id,
        "dimension": dimension,
        "number_ess_expected": number_ess,
        "is_cs": is_cs,
        "matrix": matrix_str,
        "result": result_without_candidates
    }
    
    return (matrix_id, dimension, number_ess, result_entry, candidates_list, computation_result)


def main():
    # Paths
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    matrices_file = script_dir / "verification" / "verification_matrices.json"
    
    # Generate timestamped output filenames in results folder
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_dir = script_dir / "results"
    results_dir.mkdir(exist_ok=True)
    results_file = results_dir / f"fracessa_verification_result_{timestamp}.json"
    candidates_file = results_dir / f"fracessa_verification_candidates_{timestamp}.csv"
    baseline_candidates_file = script_dir / "verification" / "baseline_candidates.csv"
    baseline_results_file = script_dir / "verification" / "baseline_result.json"

    # Check if matrices file exists
    if not matrices_file.exists():
        print(f"Error: {matrices_file} not found")
        sys.exit(1)

    # Check if fracessa_py module is available
    try:
        import fracessa_py
        print("FRACESSA Python bindings loaded successfully")
    except ImportError as e:
        print(f"Error: Cannot import fracessa_py module: {e}")
        print("Make sure fracessa_py.py is in the Python path")
        sys.exit(1)

    # Set the executable path explicitly
    executable_path = (project_root / "fracessa" / "build" / "fracessa").resolve()

    if not executable_path.exists():
        print(f"Error: fracessa executable not found at {executable_path}")
        sys.exit(1)

    # Load matrices
    print(f"Loading matrices from {matrices_file}")
    with open(matrices_file, 'r') as f:
        data = json.load(f)

    matrices = data.get('matrices', [])

    # Filter to only matrices marked as in_use
    original_count = len(matrices)
    matrices = [m for m in matrices if m.get('in_use', True)]
    skipped_count = original_count - len(matrices)

    print(f"Found {original_count} matrices total, skipping {skipped_count} matrices not in use")
    print(f"Processing {len(matrices)} matrices")

    # Sort matrices by ID for consistent processing order
    matrices.sort(key=lambda x: x.get('id', 0))
    print(f"Sorted matrices by ID: {[(m['id'], m['dimension']) for m in matrices[:5]]}...")

    # Process matrices sequentially
    results = []
    all_candidates = []
    successful = 0
    failed = 0
    total_matrices = len(matrices)

    print(f"\nStarting sequential processing...")
    start_time = time.time()

    for i, matrix_data in enumerate(matrices):
        result_tuple = process_matrix(matrix_data, str(executable_path))
        matrix_id, dimension, number_ess, result_entry, candidates_list, computation_result = result_tuple
        
        results.append(result_entry)
        all_candidates.extend(candidates_list)
        
        if computation_result["success"]:
            actual_ess = computation_result["ess_count"]
            timing = computation_result["timing"]
            status = f"✅ {actual_ess} ESS in {timing:.4f}s"
            if actual_ess != number_ess:
                status += f" ⚠️ (expected {number_ess})"
            print(f"[{i+1}/{total_matrices}] Matrix ID {matrix_id} (dim {dimension}): {status}")
            successful += 1
        else:
            error_msg = computation_result.get('error', 'Unknown error')
            if "Timed out" in error_msg:
                print(f"[{i+1}/{total_matrices}] Matrix ID {matrix_id} (dim {dimension}): ⏰ {error_msg}")
            else:
                print(f"[{i+1}/{total_matrices}] Matrix ID {matrix_id} (dim {dimension}): ❌ {error_msg}")
            failed += 1

    elapsed_time = time.time() - start_time
    print(f"\nProcessing completed in {elapsed_time:.2f}s")

    # Sort results by matrix ID for consistent output
    results.sort(key=lambda x: x['id'])
    all_candidates.sort(key=lambda x: (x['matrix_id'], x['candidate_id']))

    # Save JSON results (without candidates)
    print(f"\nSaving results to {results_file}")
    with open(results_file, 'w') as f:
        json.dump({
            "metadata": {
                "total_matrices": len(matrices),
                "successful": successful,
                "failed": failed,
                "timestamp": time.time(),
                "processing_time": elapsed_time,
                "num_processes": 1,
                "fracessa_settings": {
                    "include_candidates": True,
                    "enable_logging": False,
                    "exact_arithmetic": False,
                    "full_support_search": False,
                    "timeout": 1800.0
                }
            },
            "results": results
        }, f, indent=2)

    # Save CSV with all candidates
    print(f"Saving candidates to {candidates_file}")
    csv_columns = [
        "matrix_id", "candidate_id", "vector", "support", "support_size",
        "extended_support", "extended_support_size", "shift_reference",
        "is_ess", "reason_ess", "payoff", "payoff_double"
    ]
    
    with open(candidates_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=csv_columns)
        writer.writeheader()
        for candidate in all_candidates:
            # Convert vector list to comma-separated string
            row = candidate.copy()
            if 'vector' in row and isinstance(row['vector'], list):
                row['vector'] = ','.join(str(v) for v in row['vector'])
            writer.writerow(row)

    print("\nResults saved!")
    print(f"  JSON: {results_file}")
    print(f"  CSV:  {candidates_file}")
    print(f"\nSummary: {successful} successful, {failed} failed out of {len(matrices)} matrices")
    print(f"Total candidates: {len(all_candidates)}")
    if failed > 0:
        print("\nFailed matrices:")
        for result in results:
            if not result["result"]["success"]:
                print(f"  ID {result['id']}: {result['result']['error']}")

    # Compare ESS counts from JSON results
    print("\n" + "="*60)
    print("ESS COUNT VERIFICATION (from JSON results)")
    print("="*60)
    
    ess_count_comparison = compare_ess_counts(results, baseline_results_file)
    
    if ess_count_comparison["status"] == "SKIPPED":
        print(f"⚠️  {ess_count_comparison['reason']}")
    elif ess_count_comparison["status"] == "ERROR":
        print(f"❌ {ess_count_comparison['reason']}")
    else:
        print(f"Baseline file: {baseline_results_file}")
        print(f"Matrices compared: {ess_count_comparison['compared']}")
        
        if ess_count_comparison['compared'] > 0:
            mismatches = ess_count_comparison.get('mismatches', [])
            if mismatches:
                print(f"\n❌ ESS count mismatches found: {len(mismatches)}")
                print(f"{'Matrix ID':<10} {'Dimension':<10} {'Baseline ESS':<15} {'Current ESS':<15}")
                print("-" * 50)
                for mismatch in sorted(mismatches, key=lambda x: x['matrix_id']):
                    print(f"{mismatch['matrix_id']:<10} {mismatch['dimension']:<10} {mismatch['baseline_ess']:<15} {mismatch['current_ess']:<15}")
            else:
                print(f"\n✅ All ESS counts match baseline!")
    
    print("="*60)

    # Compare runtimes with baseline
    print("\n" + "="*60)
    print("RUNTIME COMPARISON")
    print("="*60)
    
    runtime_comparison = compare_runtimes(results, baseline_results_file)
    
    if runtime_comparison["status"] == "SKIPPED":
        print(f"⚠️  {runtime_comparison['reason']}")
    elif runtime_comparison["status"] == "ERROR":
        print(f"❌ {runtime_comparison['reason']}")
    else:
        print(f"Baseline file: {baseline_results_file}")
        print(f"Matrices compared: {runtime_comparison['compared']}")
        
        if runtime_comparison['compared'] > 0:
            # Sort details by matrix ID
            sorted_details = sorted(runtime_comparison['details'], key=lambda x: x['matrix_id'])
            
            print(f"\nRuntime comparison (current - baseline, in seconds):")
            print(f"{'Matrix ID':<10} {'Dimension':<10} {'Baseline (s)':<15} {'Current (s)':<15} {'Difference (s)':<15} {'Change %':<12}")
            print("-" * 92)
            
            for detail in sorted_details:
                matrix_id = detail['matrix_id']
                dimension = detail['dimension']
                baseline_timing = detail['baseline_timing']
                current_timing = detail['current_timing']
                difference = current_timing - baseline_timing
                percentage_change = detail.get('percentage_change')
                
                # Format difference with sign (3-4 digits after decimal for milliseconds)
                if difference > 0:
                    diff_str = f"+{difference:.4f}"
                elif difference < 0:
                    diff_str = f"{difference:.4f}"
                else:
                    diff_str = "0.0000"
                
                # Format percentage change with sign
                if percentage_change is None:
                    pct_str = "N/A"
                elif percentage_change == float('inf'):
                    pct_str = "∞"
                elif percentage_change > 0:
                    pct_str = f"+{percentage_change:.2f}%"
                elif percentage_change < 0:
                    pct_str = f"{percentage_change:.2f}%"
                else:
                    pct_str = "0.00%"
                
                print(f"{matrix_id:<10} {dimension:<10} {baseline_timing:<15.4f} {current_timing:<15.4f} {diff_str:<15} {pct_str:<12}")
    
    print("="*60)
    
    # Verify against baseline
    print("\n" + "="*60)
    print("BASELINE VERIFICATION")
    print("="*60)
    
    verification = verify_against_baseline(all_candidates, baseline_candidates_file)
    
    if verification["status"] == "SKIPPED":
        print(f"⚠️  {verification['reason']}")
    else:
        print(f"Baseline file: {baseline_candidates_file}")
        print(f"Total candidates produced: {len(all_candidates)}")
        print(f"Matrices verified: {verification['verified']}")
        print(f"Passed: {verification['passed']}")
        print(f"Failed: {verification['failed']}")
        
        # Show skipped matrices (in baseline but not in testset)
        if verification["skipped"]:
            print(f"\nℹ️  Skipped (not in testset): {len(verification['skipped'])} matrices")
            if len(verification['skipped']) > 0 and len(verification['skipped']) <= 10:
                print(f"  Matrix IDs: {verification['skipped']}")
            elif len(verification['skipped']) > 10:
                print(f"  Matrix IDs: {verification['skipped'][:10]}... (and {len(verification['skipped']) - 10} more)")
            
            # If all matrices are skipped, it might mean no candidates were produced
            if len(verification['skipped']) == len([m for m in data.get('matrices', []) if m.get('in_use', True)]):
                print(f"\n⚠️  WARNING: All matrices were skipped. This usually means:")
                print(f"    - No candidates were produced (all runs failed or returned 0 candidates)")
                print(f"    - Or matrix_id field is missing/incorrect in produced candidates")
        
        if verification["status"] == "PASS":
            print(f"\n✅ VERIFICATION PASSED - All candidates match baseline!")
        else:
            print(f"\n❌ VERIFICATION FAILED:")
            
            # Show errors (matrices not in baseline)
            if verification["errors"]:
                print("\n  Matrices not in baseline (ERROR):")
                for error in verification["errors"]:
                    print(f"    Matrix ID {error['matrix_id']}: {error['error']}")
            
            # Show candidate count mismatches
            if verification["candidate_count_mismatches"]:
                print("\n  Candidate count mismatches:")
                for mismatch in verification["candidate_count_mismatches"]:
                    print(f"    Matrix ID {mismatch['matrix_id']}: "
                          f"produced {mismatch['produced_count']}, baseline {mismatch['baseline_count']}")
            
            # Show ESS count mismatches
            if verification["ess_count_mismatches"]:
                print("\n  ESS count mismatches:")
                for mismatch in verification["ess_count_mismatches"]:
                    print(f"    Matrix ID {mismatch['matrix_id']}: "
                          f"produced {mismatch['produced_ess_count']} ESS, baseline {mismatch['baseline_ess_count']} ESS")
            
            # Show candidate set mismatches
            if verification["mismatches"]:
                print("\n  Candidate set mismatches (excluding reason_ess/stability):")
                for mismatch in verification["mismatches"]:
                    print(f"    Matrix ID {mismatch['matrix_id']}: "
                          f"produced {mismatch['produced_count']}, baseline {mismatch['baseline_count']} "
                          f"(+{mismatch['extra_count']} extra, -{mismatch['missing_count']} missing)")
    
    print("="*60)


if __name__ == "__main__":
    main()
