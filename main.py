from lagrange_method.main import LagrangeProcessor

from typing import Optional
import sympy as sp


class NewtonProcessor:
    def __init__(self, func: callable, epsilon: float, a: Optional[float] = None, b: Optional[float] = None,
                 interval_find_step: float = 1, interval_max_iterations: int = 100000,
                 interval_start_iterations: Optional[int] = None, round_count: int = 20):

        """
        :param func: The function whose root needs to be found
        :param a: Start number changing signs on an interval (if known)
        :param b: End number changing signs on an interval (if known)
        :param epsilon: Calculation accuracy float value, e.g. 10**-4 (The smaller the more accurate, but at the same time there are more iterations to solve)
        :param interval_find_step: Step value for searching intervals a, b. (it is advisable to use integers for more “human” calculations)
        :param interval_max_iterations: Maximum number of the end of iterations for searching intervals a, b
        :param interval_start_iterations: Minimum number of the start of iterations for searching intervals a, b (by default = -1interval_max_iterations)
        """

        self.x = sp.symbols('x')
        self.equations: list[str] = []

        self.function: callable = func
        self.a: Optional[float] = a
        self.b: Optional[float] = b
        self.epsilon: float = epsilon

        self.step: float = interval_find_step
        self.max_iters: int = interval_max_iterations
        self.start_iters: int = interval_start_iterations if interval_start_iterations is not None else -interval_max_iterations

        self.solutions: list[float] = []
        self.intervals: list = []

        self.working = True

        self.round_count = round_count

        self.processor = LagrangeProcessor(self.function)

    @staticmethod
    def _get_signs(integers: list[int | float]) -> str:
        return f'({", ".join(["-" if sign < 0 else "+" if sign != 0 else "" for sign in integers])})'

    def find_intervals(self) -> list[tuple[int | float, int | float]]:
        """
        :return: List of intervals where the function value changes its sign
        """
        intervals = []
        start_iters = self.start_iters

        while start_iters < self.max_iters:
            if self.function(start_iters) * self.function(start_iters + self.step) < 0:
                a = start_iters
                b = start_iters + self.step
                intervals.append((a, b))
            start_iters += self.step
        self.intervals = intervals
        return intervals

    def generate_equations(self, label: str, var: int | float, is_initial_sample: bool, c):
        """
        :param label: Label symbol in equations
        :param var: Variable of the label in equations
        :param is_initial_sample: Is initial sample?
        :param c: Result of equations
        :return: Equations
        """
        eq = f"""
    Первоначальная выборка приближения по формуле c = f(x) * f''(x) > 0.
    
    x = {label} = {var}
    f({var}) * f''({var}) = {self.function(var) * self.processor.derivative(var)} > 0
    
    Первое приближение:
    """ if is_initial_sample else ""
        eq += f"""
             f({label})               {self.function(var)}
    c = {label} - ———————— = {var} - ———————————— = {c}
             f'({label})              {self.processor.derivative(var)}
                """
        return eq

    def find_solutions(self, extra: bool = False):
        """
        :param extra: Parameter indicating whether to solve the equation using the instant method without intermediate steps
        :return: List of possible solutions
        """

        self.equations = []

        self.intervals = [(self.a, self.b)] if len(self.intervals) == 0 else self.intervals

        if extra:
            self.solutions = sp.solve(self.function(self.x), self.x)

            return self.solutions

        if self.a is None or self.b is None:
            self.intervals = self.find_intervals()

        for a, b in self.intervals:
            equation = []
            self.a = a
            self.b = b

            eq = ""

            self.working = True

            equation.append(f"\n\tПервоначально найденный интервал смены знаков: [{self.a}, {self.b}].\n")

            if self.function(self.a) * self.function(self.b) >= 0:
                raise ValueError("Функция должна иметь разные знаки на концах интервала [a, b].")

            while abs(self.b - self.a) >= self.epsilon and self.working:
                statuses = [["\033[9m", "\033[0m"], ["\033[9m", "\033[0m"]]

                if self.function(self.a) * self.processor.double_derivative(self.a) > 0 or self.a != a:
                    c = self.a - (self.function(self.a) / self.processor.derivative(self.a))
                    eq = self.generate_equations("a", self.a, len(equation) == 1, c)
                elif self.b != b:
                    c = self.b - (self.function(self.b) / self.processor.derivative(self.b))
                    eq = self.generate_equations("b", self.b, len(equation) == 1, c)

                old_a, old_b = float(self.a), float(self.b)

                if self.function(c) == 0.0:
                    equation.append(f"Ответ x = {c}. Так как f(c) = 0. (Найден точный корень)")

                    self.solutions.append(c)
                    self.equations.append(equation)
                    self.working = False

                elif self.function(c) * self.function(self.a) < 0:
                    statuses[1] = ["", ""]
                    self.b = c
                else:
                    statuses[0] = ["", ""]
                    self.a = c

                if self.working:
                    eq += f"\n\t\t{self._get_signs([self.function(old_a), self.function(c)])} {statuses[1][0]}{[old_a, c]}{statuses[1][1]} " \
                          f"| {self._get_signs([self.function(c), self.function(old_b)])} {statuses[0][0]}{[c, old_b]}{statuses[0][1]}\n"
                    equation.append(eq)

            solution = (self.a + self.b) / 2

            if self.working:
                equation.append(f"Ответ x = (a + b) / 2 = {solution}. Так как b - a < epsilon.")
                self.solutions.append(solution)
                self.equations.append(equation)
        return self.solutions

    def print_last_equations(self) -> None:
        print(f"Всего решений: {len(self.equations)}")
        for solution_index, solution in enumerate(self.equations):
            print(f"Решение {solution_index + 1}:")
            for step_index, step in enumerate(solution):
                print(f"\t{step_index + 1}. {step}")
