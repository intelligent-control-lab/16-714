# Symbolic computation in the HW2 scaffolds

The unicycle dynamics contain products of speed with trigonometric functions
of heading. Their derivatives depend on the state, so the expansion point
matters. Symbolic differentiation is useful for obtaining these expressions
before numerical evaluation. It does not replace the control derivations
required by the assignment.

## Supplied Python workflow

Question 1.4 uses SymPy in the provided `symbolic_model` and
`evaluated_jacobians` helpers:

1. `sympy.symbols(..., real=True)` declares symbolic state and control
   variables. The sampling-time symbol is also declared positive.
2. `sympy.Matrix(...)` represents the Euler dynamics as a column vector.
3. `Matrix.jacobian(...)` differentiates that vector with respect to the
   state or control vector; `sympy.simplify(...)` simplifies the expressions.
4. `subs(...)` inserts the prescribed numerical parameters and reference
   state. Symbol names and array order must stay consistent.
5. `numpy.asarray(..., dtype=float)` converts a fully evaluated expression
   into a floating-point array for the numerical backward and forward passes.

The symbolic setup is complete. Implement the mathematical TODOs using its
outputs; you do not need to modify the SymPy plumbing.

## Numerical linear systems

The first-order co-state approximation reduces the coupled equations to
small matrix systems. Use `numpy.linalg.solve(matrix, right_hand_side)` for
these systems after deriving their coefficients. This avoids explicitly
forming a numerical matrix inverse. The scaffold does not require repeated
calls to a general-purpose symbolic equation solver at every time step.

SymPy matrices and NumPy arrays have different multiplication conventions:
`*` multiplies SymPy matrices, whereas `@` performs NumPy matrix
multiplication. NumPy `*` is elementwise multiplication. Convert to numeric
arrays only after substituting every remaining symbol; otherwise the result
may contain symbolic objects instead of floating-point values.

## Questions 1.5–1.8

Question 1.5 asks you to implement the general affine linearization and LQR
steps numerically; its released scaffold does not call SymPy. Questions
1.6–1.8 also use NumPy. Question 1.7 supplies `state_jacobian` and
`delta_model_terms`, and Question 1.8 reuses them through its sibling import.
The current nominal trajectory is the input to each local-model evaluation;
no additional symbolic solver coding is needed for these questions.

The provided functions document input and output shapes. Preserve those
interfaces, the distinction between a state trajectory and a control
trajectory, and the terminal state/co-state sample. The import loader,
numerical guards, plotting, and file I/O are also provided; focus your work
on the mathematical blocks listed in [the HW2 README](README.md).
