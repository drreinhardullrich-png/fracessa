#!/usr/bin/env python3
"""
Script to run all matrices from matrices.json using the FRACESSA Python bindings
and save results to a JSON file.
"""

import json
import sys
import time
from pathlib import Path
from fracessa_py import Fracessa, Matrix, FracessaError


def run_fracessa_executable(dimension, matrix_str, executable_path, is_cs=False, include_candidates=True, enable_logging=True):
    """
    Run FRACESSA for a given matrix using the Python bindings.

    Args:
        dimension: Matrix dimension
        matrix_str: Matrix as comma-separated string
        executable_path: Path to the fracessa executable (used for Fracessa class)
        is_cs: Whether this is a circular symmetric matrix
        include_candidates: Whether to include candidate details
        enable_logging: Whether to enable logging

    Returns:
        dict: Result containing ESS count, candidates, timing, etc.
    """
    try:
        # Initialize FRACESSA interface
        fracessa = Fracessa(str(executable_path))

        # Create matrix object
        matrix = Matrix(matrix_str, dimension, is_circular=is_cs)

        # Run ESS computation
        result = fracessa.compute_ess(
            matrix=matrix,
            include_candidates=include_candidates,
            enable_logging=enable_logging,
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

            return {
                "success": True,
                "ess_count": result.ess_count,
                "timing": result.computation_time,
                "candidates": candidates_data
            }
        else:
            return {
                "success": False,
                "error": result.error,
                "timing": result.computation_time
            }

    except FracessaError as e:
        return {
            "success": False,
            "error": f"FRACESSA Error: {str(e)}",
            "timing": 0.0
        }
    except Exception as e:
        return {
            "success": False,
            "error": f"Exception: {str(e)}",
            "timing": 0.0
        }

def main():
    # Paths
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    matrices_file = script_dir / "test" / "verification_matrices.json"
    results_file = script_dir / "matrices_results.json"

    # Debug paths (optional)
    # print(f"Script directory: {script_dir}")
    # print(f"Project root: {project_root}")
    # print(f"Matrices file: {matrices_file}")
    # print(f"Results file: {results_file}")

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

    # Use sequential processing for reliability
    num_processes = 1
    print(f"Using sequential processing (1 process)")

    # Process matrices sequentially
    results = []
    successful = 0
    failed = 0

    print("Starting sequential processing...")
    total_matrices = len(matrices)
    for i, matrix in enumerate(matrices):
        matrix_id = matrix.get('id')
        dimension = matrix.get('dimension')
        number_ess = matrix.get('number_ess')
        is_cs = matrix.get('is_cs')
        matrix_str = matrix.get('matrix')
        print(f"Starting matrix {i+1}/{total_matrices}: ID {matrix_id} (dim {dimension})")
        print("Processing...", end=" ", flush=True)

        # Run fracessa with logging enabled
        result = run_fracessa_executable(
            dimension=dimension,
            matrix_str=matrix_str,
            executable_path=str(executable_path),
            is_cs=is_cs,
            include_candidates=True,
            enable_logging=True
        )

        # Build result entry
        result_entry = {
            "id": matrix_id,
            "dimension": dimension,
            "number_ess_expected": number_ess,
            "is_cs": is_cs,
            "matrix": matrix_str,
            "result": result
        }

        results.append(result_entry)

        if result["success"]:
            actual_ess = result["ess_count"]
            expected_ess = number_ess
            timing = result["timing"]
            status = f"✅ {actual_ess} ESS in {timing:.2f}s"
            if actual_ess != expected_ess:
                status += f" ⚠️ (expected {expected_ess})"
            print(f"✓ Matrix {i+1}/{total_matrices} (ID {matrix_id}, dim {dimension}): {status}")
            successful += 1
        else:
            error_msg = result['error']
            if "Timed out" in error_msg:
                print(f"⏰ Matrix {i+1}/{total_matrices} (ID {matrix_id}, dim {dimension}): {error_msg}")
            else:
                print(f"❌ Matrix {i+1}/{total_matrices} (ID {matrix_id}, dim {dimension}): Failed: {error_msg}")
            failed += 1

    # Add newline after progress display
    print()

    # Save results
    print(f"\nSaving results to {results_file}")
    with open(results_file, 'w') as f:
        json.dump({
            "metadata": {
                "total_matrices": len(matrices),
                "successful": successful,
                "failed": failed,
                "timestamp": time.time()
            },
            "results": results
        }, f, indent=2)

    print("\nResults saved!")
    print(f"Summary: {successful} successful, {failed} failed out of {len(matrices)} matrices")
    if failed > 0:
        print("\nFailed matrices:")
        for result in results:
            if not result["result"]["success"]:
                print(f"  ID {result['id']}: {result['result']['error']}")

if __name__ == "__main__":
    main()

