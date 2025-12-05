#!/usr/bin/env python3
"""
FRACESSA Python Wrapper Module

This module provides a Python API for the FRACESSA executable,
allowing you to compute evolutionary stable strategies (ESS) from Python code.
"""

import subprocess
import json
import os
from pathlib import Path
from typing import List, Dict, Optional, Union, Tuple
import time


class FracessaError(Exception):
    """Exception raised for FRACESSA-related errors."""
    pass


class Matrix:
    """
    Represents a payoff matrix for evolutionary game theory.

    Args:
        matrix_str: String representation of the matrix (e.g., "0,1,1,0" for 2x2)
        dimension: Matrix dimension (will be auto-detected if not provided)
        is_circular: Whether this is a circular symmetric matrix
    """

    def __init__(self, matrix_str: str, dimension: Optional[int] = None, is_circular: bool = False):
        self.matrix_str = matrix_str
        self.dimension = dimension
        self.is_circular = is_circular

        # Validate or auto-detect dimension
        elements = matrix_str.split(',')
        n = len(elements)
        
        if self.dimension is None:
            # Auto-detect dimension
            # Try circular symmetric first (floor(dimension/2) elements)
            # For n elements, dimension could be 2*n (even) or 2*n+1 (odd)
            # Try symmetric matrix (upper triangular: dimension*(dimension+1)/2 elements)
            # Solve: n = d*(d+1)/2 => d² + d - 2n = 0 => d = (-1 + sqrt(1+8n))/2
            import math
            symmetric_dim = int((-1 + math.sqrt(1 + 8*n)) / 2)
            if symmetric_dim * (symmetric_dim + 1) // 2 == n:
                self.dimension = symmetric_dim
                self.is_circular = False
            elif self.is_circular:
                # For circular symmetric matrices, we need dimension to be provided
                # because n elements could correspond to dimension 2*n (even) or 2*n+1 (odd)
                raise ValueError(f"Cannot auto-detect dimension for circular symmetric matrix with {n} elements. Please provide dimension.")
            else:
                raise ValueError(f"Cannot auto-detect dimension for {n} elements")
        else:
            # Dimension provided - validate number of elements
            if self.is_circular:
                # Circular symmetric: should have floor(dimension/2) elements
                expected_elements = self.dimension // 2
                if n != expected_elements:
                    raise ValueError(f"Circular symmetric matrix with dimension {self.dimension} should have {expected_elements} elements, but got {n}")
            else:
                # Symmetric matrix: should have dimension*(dimension+1)/2 elements (upper triangular)
                expected_elements = self.dimension * (self.dimension + 1) // 2
                if n != expected_elements:
                    raise ValueError(f"Symmetric matrix with dimension {self.dimension} should have {expected_elements} elements (upper triangular), but got {n}")

    def to_cli_string(self) -> str:
        """Convert to command-line format for fracessa executable."""
        return f"{self.dimension}#{self.matrix_str}"

    def __str__(self):
        return f"Matrix({self.dimension}x{self.dimension}, circular={self.is_circular})"

    def __repr__(self):
        return f"Matrix('{self.matrix_str}', dimension={self.dimension}, is_circular={self.is_circular})"


class Candidate:
    """
    Represents an ESS candidate found by FRACESSA.

    Attributes:
        candidate_id: Unique identifier for this candidate
        vector: Strategy vector (list of rational numbers as strings)
        support: Bitmask representing the support set
        support_size: Size of the support set
        extended_support: Extended support bitmask
        extended_support_size: Size of extended support
        shift_reference: Shift reference value
        is_ess: Whether this is confirmed to be an ESS
        reason_ess: Reason why this is/isn't an ESS
        payoff: Payoff value (rational)
        payoff_double: Payoff value (floating point)
    """

    def __init__(self, data: Dict):
        self.candidate_id = data.get('candidate_id', 0)
        self.vector = data.get('vector', [])
        self.support = data.get('support', 0)
        self.support_size = data.get('support_size', 0)
        self.extended_support = data.get('extended_support', 0)
        self.extended_support_size = data.get('extended_support_size', 0)
        self.shift_reference = data.get('shift_reference', 0)
        self.is_ess = data.get('is_ess', False)
        self.reason_ess = data.get('reason_ess', '')
        self.payoff = data.get('payoff', '0')
        self.payoff_double = data.get('payoff_double', 0.0)

    def __str__(self):
        return f"{'ESS' if self.is_ess else 'Candidate'} {self.candidate_id}: payoff={self.payoff_double:.6f}"

    def __repr__(self):
        return f"{'ESS' if self.is_ess else 'Candidate'}(id={self.candidate_id}, payoff={self.payoff_double:.6f}, is_ess={self.is_ess})"


class ESSResult:
    """
    Result of an ESS computation.

    Attributes:
        ess_count: Number of evolutionary stable strategies found
        candidates: List of Candidate objects
        computation_time: Time taken for computation (seconds)
        success: Whether the computation succeeded
        error: Error message if computation failed
    """

    def __init__(self, ess_count: int = 0, candidates: List[Candidate] = None,
                 computation_time: float = 0.0, success: bool = True, error: str = ""):
        self.ess_count = ess_count
        self.candidates = candidates or []
        self.computation_time = computation_time
        self.success = success
        self.error = error

    def __str__(self):
        if self.success:
            return f"ESS Result: {self.ess_count} ESS found in {self.computation_time:.2f}s"
        else:
            return f"ESS Result: Failed - {self.error}"

    def __repr__(self):
        return f"ESSResult(ess_count={self.ess_count}, candidates={len(self.candidates)}, success={self.success})"


class Fracessa:
    """
    Main FRACESSA interface for computing evolutionary stable strategies.

    Args:
        executable_path: Path to the fracessa executable (auto-detected if not provided)
    """

    def __init__(self, executable_path: Optional[str] = None):
        if executable_path:
            self.executable_path = Path(executable_path)
        else:
            # Auto-detect executable path
            script_dir = Path(__file__).parent.parent
            possible_paths = [
                script_dir / "fracessa" / "build" / "fracessa",
                script_dir / "build" / "fracessa",
                Path("./fracessa/build/fracessa"),
                Path("./build/fracessa")
            ]

            self.executable_path = None
            for path in possible_paths:
                if path.exists() and path.is_file():
                    self.executable_path = path
                    break

            if self.executable_path is None:
                raise FracessaError("Could not find fracessa executable. Please specify executable_path.")

        if not self.executable_path.exists():
            raise FracessaError(f"Fracessa executable not found at {self.executable_path}")

        # Make executable
        os.chmod(self.executable_path, 0o755)

    def compute_ess(self, matrix: Union[Matrix, str],
                   include_candidates: bool = True,
                   exact_arithmetic: bool = False,
                   full_support_search: bool = False,
                   enable_logging: bool = False,
                   timeout: float = 0.0) -> ESSResult:
        """
        Compute evolutionary stable strategies for a given payoff matrix.

        Args:
            matrix: Matrix object or CLI string (e.g., "2#0,1,1,0")
            include_candidates: Whether to return detailed candidate information
            exact_arithmetic: Use exact rational arithmetic (slower but precise)
            full_support_search: Search full support directly after size 1
            enable_logging: Enable detailed logging to fracessa.log
            timeout: Maximum computation time in seconds (0.0 means no timeout)

        Returns:
            ESSResult object containing the computation results

        Raises:
            FracessaError: If computation fails
        """

        # Convert matrix to CLI string
        if isinstance(matrix, Matrix):
            cli_string = matrix.to_cli_string()
        elif isinstance(matrix, str):
            cli_string = matrix
        else:
            raise FracessaError("matrix must be a Matrix object or CLI string")

        # Build command
        cmd = [str(self.executable_path)]
        if include_candidates:
            cmd.append("-c")
        if exact_arithmetic:
            cmd.append("-e")
        if full_support_search:
            cmd.append("-f")
        if enable_logging:
            cmd.append("-l")
        cmd.append(cli_string)

        try:
            # Run command with timeout (0.0 means no timeout)
            start_time = time.time()
            subprocess_timeout = None if timeout == 0.0 else timeout
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=subprocess_timeout
            )
            end_time = time.time()
            computation_time = end_time - start_time

            if result.returncode != 0:
                return ESSResult(
                    success=False,
                    error=f"Command failed with return code {result.returncode}: {result.stderr}",
                    computation_time=computation_time
                )

            # Parse output
            lines = result.stdout.strip().split('\n')
            if not lines or not lines[0].strip():
                return ESSResult(
                    success=False,
                    error="No output from executable",
                    computation_time=computation_time
                )

            try:
                ess_count = int(lines[0].strip())
            except ValueError:
                return ESSResult(
                    success=False,
                    error=f"Could not parse ESS count from: {lines[0]}",
                    computation_time=computation_time
                )

            candidates = []
            if include_candidates and len(lines) > 1:
                if len(lines) > 1:
                    header = lines[1]  # Header line
                    candidate_lines = lines[2:]  # Candidate data lines

                    for line in candidate_lines:
                        if line.strip():
                            parts = line.split(';')
                            if len(parts) >= 11:
                                candidate_data = {
                                    'candidate_id': int(parts[0]) if parts[0] else 0,
                                    'vector': parts[1].split(',') if parts[1] else [],
                                    'support': int(parts[2]) if parts[2] else 0,
                                    'support_size': int(parts[3]) if parts[3] else 0,
                                    'extended_support': int(parts[4]) if parts[4] else 0,
                                    'extended_support_size': int(parts[5]) if parts[5] else 0,
                                    'shift_reference': int(parts[6]) if parts[6] else 0,
                                    'is_ess': parts[7].lower() == '1' if parts[7] else False,
                                    'reason_ess': parts[8] if len(parts) > 8 else '',
                                    'payoff': parts[9] if len(parts) > 9 else '0',
                                    'payoff_double': float(parts[10]) if len(parts) > 10 and parts[10] else 0.0
                                }
                                candidates.append(Candidate(candidate_data))

            return ESSResult(
                ess_count=ess_count,
                candidates=candidates,
                computation_time=computation_time,
                success=True
            )

        except subprocess.TimeoutExpired:
            return ESSResult(
                success=False,
                error=f"Computation timed out after {timeout} seconds",
                computation_time=timeout
            )
        except Exception as e:
            return ESSResult(
                success=False,
                error=f"Exception: {str(e)}",
                computation_time=0.0
            )

    def __repr__(self):
        return f"Fracessa(executable='{self.executable_path}')"
