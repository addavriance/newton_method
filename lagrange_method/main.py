import inspect


class LagrangeProcessor:

    """
    Class that finds the solutions of the derivative and double derivative of a given equation f with a specified accuracy h.
     |
    To find the derivative and double derivative, use the class methods .derivative(x) and .double_derivative(x) respectively.
     |
    After running one of the derivative finding methods, the result will be returned. Additionally, after running and finding the result, the solution steps will be added to the class list .equations.
    """

    def __init__(self, function: callable, round_count: int = 10, h: float = 0.00001):
        """
        :param function: Function that will be processed
        :param round_count: Number of decimal places when rounding
        :param h: Fractional number accuracy of calculations.
        The lower the value, the better, but if you set the values too low(less that 10^-4), inaccuracy will increase or result will be incorrect.
        """

        self.function: callable = function
        self.round_count = round_count
        self.h: float = h

        self.equations: list[str] = []
        self.solutions: list[int | float] = []

    def get_str_func(self) -> str:
        source_code = inspect.getsource(self.function).split("return ")[-1].split("#")[0].strip()

        return source_code

    def derivative(self, x: int | float) -> int | float:
        """
        :param x: Value with which derivative will be calculated.
        :return: The result of finding the derivative of a given function in initializer with a specified argument X.
        """

        xNp1 = x + self.h
        xNm1 = x - self.h
        dX = 2 * self.h

        fNp1 = round(self.function(xNp1), self.round_count)
        fNm1 = round(self.function(xNm1), self.round_count)

        derivative = round((fNp1-fNm1)/dX, self.round_count)

        equation = f"""
    Решение производной функции {self.get_str_func()} со значением xₙ = {x}.
        
        xₙ₊₁ = xₙ + h = {x} + {self.h} = {xNp1}
        xₙ₋₁ = xₙ - h = {x} - {self.h} = {xNm1}
        Δx = 2 * h = 2 * {self.h} = {dX}
        
        Fₙ₊₁ = f(xₙ₊₁) = f({xNp1}) = {fNp1}
        Fₙ₋₁ = f(xₙ₋₁) = f({xNm1}) = {fNm1}
        
                Fₙ₊₁ - Fₙ₋₁      {fNp1} - {fNm1}
        F' = ———————————————— = ———————————————— = {derivative}
                    Δx                 {dX}
        """

        self.equations.append(equation)

        self.solutions.append(derivative)

        return derivative

    def double_derivative(self, x: int | float) -> int | float:
        """
        :param x: Value with which double derivative will be calculated.
        :return: The result of finding the double derivative of a given function in initializer with a specified argument X.
        """

        xNp2 = x + 2 * self.h
        xNm2 = x - 2 * self.h
        dX = 2 * self.h

        fN = round(self.function(x), self.round_count)
        fNp2 = round(self.function(xNp2), self.round_count)
        fNm2 = round(self.function(xNm2), self.round_count)

        double_derivative = round(((fNp2 - 2*fN + fNm2) / dX**2), self.round_count)

        equation = f"""
    Решение двойной производной функции {self.get_str_func()} со значением xₙ = {x}.
        
        xₙ₊₂ = xₙ + 2 * h = {x} + 2 * {self.h} = {xNp2}
        xₙ₋₂ = xₙ - 2 * h = {x} - 2 * {self.h} = {xNm2}
        Δx = 2 * h = 2 * {self.h} = {dX}

        Fₙ₊₂ = f(xₙ₊₂) = f({xNp2}) = {fNp2}
        Fₙ₋₂ = f(xₙ₋₂) = f({xNm2}) = {fNm2}
        Fₙ = f(xₙ) = {fN}

              Fₙ₊₂ - 2*Fₙ + Fₙ₋₂  {fNp2} - 2*{fN} + {fNm2}  
        F'' = ———————————————— = ———————————————— = {double_derivative}
                    (Δx)^2              {dX**2}
        """

        self.equations.append(equation)

        self.solutions.append(double_derivative)

        return double_derivative
