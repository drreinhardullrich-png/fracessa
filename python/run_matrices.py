#!/usr/bin/env python3
"""
Script to run all matrices from matrices.json using the FRACESSA Python bindings
and save results to a JSON file. Uses multiprocessing for parallel execution.
"""

import json
import sys
import time
import csv
import multiprocessing as mp
from datetime import datetime
from pathlib import Path
from fracessa_py import Fracessa, Matrix, FracessaError


# Global variable for executable path (needed for multiprocessing)
_executable_path = None


def init_worker(executable_path):
    """Initialize worker process with executable path."""
    global _executable_path
    _executable_path = executable_path


def candidate_to_comparison_key(candidate):
    """
    Convert a candidate dict to a hashable tuple for comparison.
    Excludes candidate_id since it may differ between runs.
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
        str(candidate.get('reason_ess', '')),
        str(candidate.get('payoff', '')),
        str(candidate.get('payoff_double', ''))
    )


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
    
    for candidate in all_candidates:
        matrix_id = candidate.get('matrix_id')
        key = candidate_to_comparison_key(candidate)
        produced_by_matrix[matrix_id].add(key)
    
    for candidate in baseline_candidates:
        matrix_id = int(candidate.get('matrix_id'))
        key = candidate_to_comparison_key(candidate)
        baseline_by_matrix[matrix_id].add(key)
    
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
        
        if produced_set == baseline_set:
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
        "mismatches": mismatches
    }


def process_matrix(matrix_data):
    """
    Worker function to process a single matrix.
    
    Args:
        matrix_data: dict with matrix information
        
    Returns:
        tuple: (matrix_id, result_entry, candidates_list)
    """
    global _executable_path
    
    matrix_id = matrix_data['id']
    dimension = matrix_data['dimension']
    number_ess = matrix_data['number_ess']
    is_cs = matrix_data['is_cs']
    matrix_str = matrix_data['matrix']
    
    try:
        # Initialize FRACESSA interface
        fracessa = Fracessa(str(_executable_path))
        
        # Create matrix object
        matrix = Matrix(matrix_str, dimension, is_circular=is_cs)
        
        # Run ESS computation
        result = fracessa.compute_ess(
            matrix=matrix,
            include_candidates=True,
            enable_logging=True,
            timeout=1800.0  # 30 minutes timeout
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
    matrices_file = script_dir / "test" / "verification_matrices.json"
    
    # Generate timestamped output filenames in results folder
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    results_dir = script_dir / "results"
    results_dir.mkdir(exist_ok=True)
    results_file = results_dir / f"fracessa_verification_result_{timestamp}.json"
    candidates_file = results_dir / f"fracessa_verification_candidates_{timestamp}.csv"
    baseline_file = script_dir / "test" / "fracessa_verification_candidates_baseline.csv"

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

    # Sort matrices by dimension in descending order (largest first)
    matrices.sort(key=lambda x: x.get('dimension', 0), reverse=True)
    print(f"Sorted matrices by dimension (largest first): {[(m['id'], m['dimension']) for m in matrices[:5]]}...")

    # Determine number of processes (all cores except one)
    num_cores = mp.cpu_count()
    num_processes = max(1, num_cores - 1)
    print(f"Using {num_processes} processes ({num_cores} cores available)")

    # Process matrices in parallel
    results = []
    all_candidates = []
    successful = 0
    failed = 0
    total_matrices = len(matrices)

    print(f"\nStarting parallel processing with {num_processes} workers...")
    start_time = time.time()

    with mp.Pool(processes=num_processes, initializer=init_worker, initargs=(str(executable_path),)) as pool:
        # Use imap_unordered for results as they complete
        for i, result_tuple in enumerate(pool.imap_unordered(process_matrix, matrices)):
            matrix_id, dimension, number_ess, result_entry, candidates_list, computation_result = result_tuple
            
            results.append(result_entry)
            all_candidates.extend(candidates_list)
            
            if computation_result["success"]:
                actual_ess = computation_result["ess_count"]
                timing = computation_result["timing"]
                status = f"✅ {actual_ess} ESS in {timing:.2f}s"
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
                "num_processes": num_processes,
                "fracessa_settings": {
                    "include_candidates": True,
                    "enable_logging": True,
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

    # Verify against baseline
    print("\n" + "="*60)
    print("BASELINE VERIFICATION")
    print("="*60)
    
    verification = verify_against_baseline(all_candidates, baseline_file)
    
    if verification["status"] == "SKIPPED":
        print(f"⚠️  {verification['reason']}")
    else:
        print(f"Baseline file: {baseline_file}")
        print(f"Matrices verified: {verification['verified']}")
        print(f"Passed: {verification['passed']}")
        print(f"Failed: {verification['failed']}")
        
        # Show skipped matrices (in baseline but not in testset)
        if verification["skipped"]:
            print(f"\nℹ️  Skipped (not in testset): {verification['skipped']}")
        
        if verification["status"] == "PASS":
            print(f"\n✅ VERIFICATION PASSED - All candidates match baseline!")
        else:
            print(f"\n❌ VERIFICATION FAILED:")
            
            # Show errors (matrices not in baseline)
            if verification["errors"]:
                print("\n  Matrices not in baseline (ERROR):")
                for error in verification["errors"]:
                    print(f"    Matrix ID {error['matrix_id']}: {error['error']}")
            
            # Show mismatches
            if verification["mismatches"]:
                print("\n  Candidate mismatches:")
                for mismatch in verification["mismatches"]:
                    print(f"    Matrix ID {mismatch['matrix_id']}: "
                          f"produced {mismatch['produced_count']}, baseline {mismatch['baseline_count']} "
                          f"(+{mismatch['extra_count']} extra, -{mismatch['missing_count']} missing)")
    
    print("="*60)


if __name__ == "__main__":
    main()
