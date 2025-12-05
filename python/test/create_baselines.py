#!/usr/bin/env python3
"""
Script to regenerate baseline files by running the old executable (REF_linux64.exe)
on all matrices from verification_matrices.json.

This script:
1. Loads matrices from verification_matrices.json
2. Runs the old executable (python/the_old_exe/REF_linux64.exe) with -c flag
3. Parses the output to extract ESS counts and candidates
4. Writes both baseline files:
   - fracessa_verification_result_baseline.json
   - fracessa_verification_candidates_baseline.csv
"""

import json
import csv
import subprocess
import sys
import time
import os
from pathlib import Path
from typing import List, Dict, Tuple, Optional


def full_matrix_to_upper_triangular(matrix_str: str, dimension: int) -> str:
    """
    Convert a full matrix (n*n elements in row-major order) to upper triangular format.
    
    For an n×n matrix, upper triangular includes elements where i <= j.
    In row-major order: row i starts at index i*n, and we take elements from column i to n-1.
    
    Args:
        matrix_str: Comma-separated string of n*n elements (row-major order)
        dimension: Matrix dimension n
        
    Returns:
        Comma-separated string of upper triangular elements
    """
    elements = matrix_str.split(',')
    if len(elements) != dimension * dimension:
        raise ValueError(f"Expected {dimension * dimension} elements, got {len(elements)}")
    
    upper_tri_elements = []
    for i in range(dimension):
        # Row i: take elements from column i to n-1
        row_start = i * dimension
        for j in range(i, dimension):
            upper_tri_elements.append(elements[row_start + j])
    
    return ','.join(upper_tri_elements)


def find_old_executable(script_dir: Path) -> Path:
    """Find the old executable path."""
    # Try relative to script directory (python/test/)
    project_root = script_dir.parent.parent
    old_exe_path = project_root / "python" / "the_old_exe" / "REF_linux64.exe"
    
    if old_exe_path.exists():
        return old_exe_path
    
    # Try absolute path
    old_exe_path = Path("python/the_old_exe/REF_linux64.exe")
    if old_exe_path.exists():
        return old_exe_path.resolve()
    
    raise FileNotFoundError(f"Old executable not found at {old_exe_path}")


def run_old_executable(exe_path: Path, matrix_cli_string: str, timeout: float = 1800.0) -> Tuple[bool, int, List[Dict], float, Optional[str]]:
    """
    Run the old executable and parse output.
    
    Returns:
        (success, ess_count, candidates, timing, error_message)
    """
    cmd = [str(exe_path), "-c", matrix_cli_string]
    
    try:
        start_time = time.time()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout
        )
        end_time = time.time()
        timing = end_time - start_time
        
        if result.returncode != 0:
            return (False, 0, [], timing, f"Command failed with return code {result.returncode}: {result.stderr}")
        
        # Parse output
        lines = result.stdout.strip().split('\n')
        if not lines or not lines[0].strip():
            return (False, 0, [], timing, "No output from executable")
        
        try:
            ess_count = int(lines[0].strip())
        except ValueError:
            return (False, 0, [], timing, f"Could not parse ESS count from: {lines[0]}")
        
        # Parse candidates
        candidates = []
        if len(lines) > 1:
            # Second line is header, skip it
            candidate_lines = lines[2:] if len(lines) > 2 else []
            
            for line in candidate_lines:
                if line.strip():
                    parts = line.split(';')
                    if len(parts) >= 11:
                        # Parse is_ess: could be "1", "True", "true", etc.
                        is_ess_str = parts[7].strip().lower() if len(parts) > 7 else "0"
                        is_ess = is_ess_str in ("1", "true", "yes")
                        
                        candidate_data = {
                            'candidate_id': int(parts[0]) if parts[0] else 0,
                            'vector': parts[1] if parts[1] else '',  # Keep as string (comma-separated)
                            'support': int(parts[2]) if parts[2] else 0,
                            'support_size': int(parts[3]) if parts[3] else 0,
                            'extended_support': int(parts[4]) if parts[4] else 0,
                            'extended_support_size': int(parts[5]) if parts[5] else 0,
                            'shift_reference': int(parts[6]) if parts[6] else 0,
                            'is_ess': is_ess,
                            'reason_ess': parts[8] if len(parts) > 8 else '',
                            'payoff': parts[9] if len(parts) > 9 else '0',
                            'payoff_double': float(parts[10]) if len(parts) > 10 and parts[10] else 0.0
                        }
                        candidates.append(candidate_data)
        
        return (True, ess_count, candidates, timing, None)
        
    except subprocess.TimeoutExpired:
        return (False, 0, [], timeout, f"Computation timed out after {timeout} seconds")
    except Exception as e:
        return (False, 0, [], 0.0, f"Exception: {str(e)}")


def process_matrix(matrix_data: Dict, old_exe_path: Path) -> Tuple[Dict, List[Dict]]:
    """
    Process a single matrix.
    
    Returns:
        (result_entry, candidates_list)
    """
    matrix_id = matrix_data['id']
    dimension = matrix_data['dimension']
    number_ess = matrix_data['number_ess']
    is_cs = matrix_data['is_cs']
    matrix_str = matrix_data['matrix']
    
    # Build CLI string: dimension#matrix_string
    # For circular symmetric: matrix_str already has n/2 elements
    # For non-circular symmetric: matrix_str already has n*n elements (full matrix in row-major order)
    cli_string = f"{dimension}#{matrix_str}"
    
    # Run old executable
    success, ess_count, candidates, timing, error = run_old_executable(old_exe_path, cli_string)
    
    if success:
        # For non-circular symmetric matrices, convert to upper triangular format for storage
        # (but keep full matrix format for CLI string when calling old executable)
        if not is_cs:
            matrix_str_for_output = full_matrix_to_upper_triangular(matrix_str, dimension)
        else:
            matrix_str_for_output = matrix_str
        
        result_entry = {
            "id": matrix_id,
            "dimension": dimension,
            "number_ess_expected": number_ess,
            "is_cs": is_cs,
            "matrix": matrix_str_for_output,
            "result": {
                "success": True,
                "ess_count": ess_count,
                "timing": timing
            }
        }
        
        # Add matrix_id to each candidate
        candidates_list = []
        for candidate in candidates:
            candidate_with_matrix_id = {"matrix_id": matrix_id}
            candidate_with_matrix_id.update(candidate)
            candidates_list.append(candidate_with_matrix_id)
        
        return (result_entry, candidates_list)
    else:
        # For non-circular symmetric matrices, convert to upper triangular format for storage
        if not is_cs:
            matrix_str_for_output = full_matrix_to_upper_triangular(matrix_str, dimension)
        else:
            matrix_str_for_output = matrix_str
        
        result_entry = {
            "id": matrix_id,
            "dimension": dimension,
            "number_ess_expected": number_ess,
            "is_cs": is_cs,
            "matrix": matrix_str_for_output,
            "result": {
                "success": False,
                "error": error,
                "timing": timing
            }
        }
        return (result_entry, [])


def main():
    """Main function."""
    script_dir = Path(__file__).resolve().parent
    matrices_file = script_dir / "verification_matrices.json"
    baseline_results_file = script_dir / "fracessa_verification_result_baseline.json"
    baseline_candidates_file = script_dir / "fracessa_verification_candidates_baseline.csv"
    
    # Find old executable
    try:
        old_exe_path = find_old_executable(script_dir)
        print(f"Using old executable: {old_exe_path}")
        
        # Make executable if needed
        os.chmod(old_exe_path, 0o755)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    # Load matrices
    if not matrices_file.exists():
        print(f"Error: {matrices_file} not found")
        sys.exit(1)
    
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
    print()
    
    # Process matrices sequentially
    results = []
    all_candidates = []
    successful = 0
    failed = 0
    total_matrices = len(matrices)
    
    start_time = time.time()
    
    for i, matrix_data in enumerate(matrices):
        matrix_id = matrix_data['id']
        dimension = matrix_data['dimension']
        
        print(f"[{i+1}/{total_matrices}] Processing Matrix ID {matrix_id} (dim {dimension})...", end=" ", flush=True)
        
        result_entry, candidates_list = process_matrix(matrix_data, old_exe_path)
        results.append(result_entry)
        all_candidates.extend(candidates_list)
        
        if result_entry["result"]["success"]:
            actual_ess = result_entry["result"]["ess_count"]
            timing = result_entry["result"]["timing"]
            expected_ess = matrix_data['number_ess']
            status = f"✅ {actual_ess} ESS in {timing:.2f}s"
            if actual_ess != expected_ess:
                status += f" ⚠️ (expected {expected_ess})"
            print(status)
            successful += 1
        else:
            error_msg = result_entry["result"].get('error', 'Unknown error')
            print(f"❌ {error_msg}")
            failed += 1
    
    elapsed_time = time.time() - start_time
    print()
    print(f"Processing completed in {elapsed_time:.2f}s")
    
    # Sort results by matrix ID for consistent output
    results.sort(key=lambda x: x['id'])
    all_candidates.sort(key=lambda x: (x['matrix_id'], x['candidate_id']))
    
    # Write JSON baseline file
    print(f"\nWriting results to {baseline_results_file}")
    with open(baseline_results_file, 'w') as f:
        json.dump({
            "metadata": {
                "total_matrices": len(matrices),
                "successful": successful,
                "failed": failed,
                "timestamp": time.time(),
                "processing_time": elapsed_time,
                "num_processes": 1,  # Sequential processing
                "fracessa_settings": {
                    "include_candidates": True,
                    "enable_logging": False,  # Old executable may not support this
                    "exact_arithmetic": False,
                    "full_support_search": False,
                    "timeout": 1800.0
                }
            },
            "results": results
        }, f, indent=2)
    
    # Write CSV baseline file
    print(f"Writing candidates to {baseline_candidates_file}")
    csv_columns = [
        "matrix_id", "candidate_id", "vector", "support", "support_size",
        "extended_support", "extended_support_size", "shift_reference",
        "is_ess", "reason_ess", "payoff", "payoff_double"
    ]
    
    with open(baseline_candidates_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=csv_columns)
        writer.writeheader()
        for candidate in all_candidates:
            # Convert is_ess boolean to string representation for CSV
            row = candidate.copy()
            row['is_ess'] = 'True' if row.get('is_ess', False) else 'False'
            writer.writerow(row)
    
    print("\n✅ Baseline files created successfully!")
    print(f"  JSON: {baseline_results_file}")
    print(f"  CSV:  {baseline_candidates_file}")
    print(f"\nSummary: {successful} successful, {failed} failed out of {len(matrices)} matrices")
    print(f"Total candidates: {len(all_candidates)}")
    
    if failed > 0:
        print("\nFailed matrices:")
        for result in results:
            if not result["result"]["success"]:
                print(f"  ID {result['id']}: {result['result']['error']}")


if __name__ == "__main__":
    main()

