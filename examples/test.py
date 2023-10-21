from main import NewtonProcessor


def f(x: int | float) -> int | float:
    return x**2-x**3-32  # replace with your function


epsilon = 10 ** -3

processor = NewtonProcessor(f, epsilon=epsilon)

processor.find_solutions()
processor.print_last_equations()
