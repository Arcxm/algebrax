from math import *
import cmath as cm

class Vec3:
    """
    A class representing a 3-dimensional vector in R^3.

    Attributes:
        x (float): x component
        y (float): y component
        z (float): z component
    """

    def __init__(self: 'Vec3', x: float, y: float, z: float) -> 'Vec3':
        """
        Initializes a new Vec3 with the given components.

        Parameters:
            x (float): x component
            y (float): y component
            z (float): z component
        
        Returns:
            Vec3: The new Vec3
        """
        self.x = x
        self.y = y
        self.z = z
    
    def __str__(self: 'Vec3') -> str:
       return f"({self.x}, {self.y}, {self.z})"
    
    def __repr__(self: 'Vec3') -> str:
       return f"Vec3({self.x}, {self.y}, {self.z})"
    
    def __eq__(self: 'Vec3', rhs) -> bool:
        if isinstance(rhs, Vec3):
            return (self.x == rhs.x and self.y == rhs.y and self.z == rhs.z)
        else:
            return NotImplemented

    def __ne__(self: 'Vec3', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq

    def __abs__(self: 'Vec3') -> float:
        return self.norm()

    def __neg__(self: 'Vec3') -> 'Vec3':
        return Vec3(-self.x, -self.y, -self.z)

    def __add__(self: 'Vec3', rhs) -> 'Vec3':
        if isinstance(rhs, Vec3):
            return Vec3(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
        else:
            return NotImplemented
    
    def __radd__(self: 'Vec3', lhs) -> 'Vec3':
        if isinstance(lhs, Vec3):
            return Vec3(lhs.x + self.x, lhs.y + self.y, lhs.z + self.z)
        else:
            return NotImplemented
    
    def __sub__(self: 'Vec3', rhs) -> 'Vec3':
        if isinstance(rhs, Vec3):
            return Vec3(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
        else:
            return NotImplemented
    
    def __rsub__(self: 'Vec3', lhs) -> 'Vec3':
        if isinstance(lhs, Vec3):
            return Vec3(lhs.x - self.x, lhs.y - self.y, lhs.z - self.z)
        else:
            return NotImplemented
    
    def __mul__(self: 'Vec3', rhs) -> 'Vec3':
        if isinstance(rhs, (int, float)):
            return Vec3(self.x * rhs, self.y * rhs, self.z * rhs)
        else:
            return NotImplemented
    
    def __rmul__(self: 'Vec3', lhs) -> 'Vec3':
        if isinstance(lhs, (int, float)):
            return Vec3(lhs * self.x, lhs * self.y, lhs * self.z)
        else:
            return NotImplemented
    
    def __truediv__(self: 'Vec3', rhs) -> 'Vec3':
        if isinstance(rhs, (int, float)):
            return Vec3(self.x / rhs, self.y / rhs, self.z / rhs)
        else:
            return NotImplemented

    def norm(self: 'Vec3') -> float:
        """
        Calculates the norm of the Vec3.
        
            sqrt(x^2 + y^2 + z^2)

        Returns:
            float: The norm of this Vec3
        """
        return sqrt(self.x**2 + self.y**2 + self.z**2)

    def normalize(self: 'Vec3') -> 'Vec3':
        """
        Normalizes the Vec3 so that its length equals 1.

        Returns:
            Vec3: A Vec3 with the same direction but length 1
        """
        n = self.norm()
        return Vec3(self.x / n, self.y / n, self.z / n)

    def dot(self: 'Vec3', other: 'Vec3') -> float:
        """
        Calculates the dot-product of this Vec3 with another Vec3.

            (x_1, y_1, z_1) * (x_2, y_2, z_2) = x_1 * x_2 + y_1 * y_2 + z_1 * z_2

        Parameters:
            other (Vec3): other Vec3
        
        Returns:
            float: The dot-product of this Vec3 with the other Vec3
        """
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self: 'Vec3', other: 'Vec3') -> 'Vec3':
        """
        Calculates the cross-product of this Vec3 with another Vec3.

            (x_1, y_1, z_1) x (x_2, y_2, z_2) = (y_1 * z_2 - z_1 * y_2, z_1 * x_2 - x_1 * z_2, x_1 * y_2 - y_1 * x_2)

        Parameters:
            other (Vec3): other Vec3
        
        Returns:
            Vec3: The cross-product of this Vec3 with the other Vec3
        """
        return Vec3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x
        )

class Dual:
    """
    A class representing a dual number.

    A dual number is a number of the form
        a + b * ϵ
    with a, b in R and ϵ^2 = 0.

    Attributes:
        a (float): real component
        b (float): dual component
    """
    
    def __init__(self: 'Dual', a: float, b: float) -> 'Dual':
        """
        Initializes a new dual number with the given components.

        Parameters:
            a (float): real component
            b (float): dual component
        
        Returns:
            Dual: The new dual number
        """
        self.a = a
        self.b = b
    
    def __str__(self: 'Dual') -> str:
        return f"({self.a}, {self.b})"
    
    def __repr__(self: 'Dual') -> str:
        return f"Dual({self.a}, {self.b})"

    def __eq__(self: 'Dual', rhs) -> bool:
        if isinstance(rhs, (int, float)):
            return (self.a == rhs and self.b == 0)
        elif isinstance(rhs, Dual):
            return (self.a == rhs.a and self.b == rhs.b)
        else:
            return NotImplemented

    def __ne__(self: 'Dual', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq

    def __abs__(self: 'Dual') -> float:
        return self.norm()
    
    def __neg__(self: 'Dual') -> 'Dual':
        return Dual(-self.a, -self.b)

    def __add__(self: 'Dual', rhs) -> 'Dual':
        if isinstance(rhs, (int, float)):
            return Dual(self.a + rhs, self.b)
        elif isinstance(rhs, Dual):
            return Dual(self.a + rhs.a, self.b + rhs.b)
        else:
            return NotImplemented

    def __radd__(self: 'Dual', lhs) -> 'Dual':
        if isinstance(lhs, (int, float)):
            return Dual(lhs + self.a, self.b)
        elif isinstance(lhs, Dual):
            return Dual(lhs.a + self.a, lhs.b + self.b)
        else:
            return NotImplemented
    
    def __sub__(self: 'Dual', rhs) -> 'Dual':
        if isinstance(rhs, (int, float)):
            return Dual(self.a - rhs, self.b)
        elif isinstance(rhs, Dual):
            return Dual(self.a - rhs.a, self.b - rhs.b)
        else:
            return NotImplemented
    
    def __rsub__(self: 'Dual', lhs) -> 'Dual':
        if isinstance(lhs, (int, float)):
            return Dual(lhs - self.a, -self.b)
        elif isinstance(lhs, Dual):
            return Dual(lhs.a - self.a, lhs.b - self.b)
        else:
            return NotImplemented
    
    def __mul__(self: 'Dual', rhs) -> 'Dual':
        if isinstance(rhs, (int, float)):
            return Dual(self.a * rhs, self.b * rhs)
        elif isinstance(rhs, Dual):
            return Dual(self.a * rhs.a, self.a * rhs.b + rhs.a * self.b)
        else:
            return NotImplemented
    
    def __rmul__(self: 'Dual', lhs) -> 'Dual':
        if isinstance(lhs, (int, float)):
            return Dual(lhs * self.a, lhs * self.b)
        elif isinstance(lhs, Dual):
            return Dual(lhs.a * self.a, lhs.a * self.b + self.a * lhs.b)
        else:
            return NotImplemented
    
    def __truediv__(self: 'Dual', rhs) -> 'Dual':
        if isinstance(rhs, (int, float)):
            return Dual(self.a / rhs, self.b / rhs)
        elif isinstance(rhs, Dual):
            return Dual(self.a / rhs.a, (self.b * rhs.a - self.a * rhs.b) / rhs.a**2)
        else:
            return NotImplemented
    
    def __rtruediv__(self: 'Dual', lhs) -> 'Dual':
        if isinstance(lhs, (int, float)):
            return Dual(lhs / self.a, (-lhs * self.b) / self.a**2)
        elif isinstance(lhs, Dual):
            return Dual(lhs.a / self.a, (lhs.b * self.a - lhs.a * self.b) / self.a**2)
        else:
            return NotImplemented

    def __pow__(self: 'Dual', exp) -> 'Dual':
        if isinstance(exp, int):
            z = Dual(1, 0)

            if exp == 0:
                return z
            
            mul = self if exp > 0 else self.inverse()

            for _ in range(abs(exp)):
                z *= mul
            return z
        elif isinstance(exp, Dual):
            if self.a <= 0:
                raise ValueError("Pow undefined: base dual number requires a > 0")
            else:
                return Dual.exp(exp * Dual.log(self))
        else:
            return NotImplemented
    
    def __rpow__(self: 'Dual', base) -> 'Dual':
        if isinstance(base, (int, float)):
            if base <= 0:
                raise ValueError("Pow undefined: base requires x > 0")
            else:
                z = Dual(base, 0)
                return z**self
        elif isinstance(base, Dual):
            return base**self
        else:
            return NotImplemented
    
    def conj(self: 'Dual') -> 'Dual':
        """
        Returns the conjugate of this dual number.

        The conjugate of a dual number
            a + b * ϵ
        is
            a - b * ϵ

        Returns:
            Dual: The conjugate of this dual number
        """
        return Dual(self.a, -self.b)

    def norm(self: 'Dual') -> float:
        """
        Calculates the norm of this dual number.
        
            a^2

        Returns:
            float: The norm of this dual number
        """
        return self.a**2

    def inverse(self: 'Dual') -> 'Dual':
        """
        Calculates the inverse of this dual number.

            1/a - b/a^2 * ϵ
        
        Returns:
            Dual: The inverse of this dual number
        """
        return Dual(1 / self.a, -self.b / self.a**2)
    
    @staticmethod
    def E() -> 'Dual':
        """
        Returns
            0 + 1 * ϵ

        Returns:
            Dual: ϵ
        """
        return Dual(0, 1)
    
    @staticmethod
    def expE(x: float) -> 'Dual':
        """
        Calculates e to the power of ϵ * x.

            e^(ϵ * x) = 1 + x * ϵ
        
        Parameters:
            x (float): coefficient of ϵ
        
        Returns:
            Dual: The dual number representing e^(ϵ * x)
        """
        return Dual(1, x)

    @staticmethod
    def exp(z: 'Dual') -> 'Dual':
        """
        Calculates e to the power of a dual number.

            e^(a + b * ϵ) = e^a * (1 + b * ϵ)
        
        Parameters:
            z (Dual): The dual number
        
        Returns:
            Dual: e raised to the power of the dual number
        """
        return exp(z.a) * Dual(1, z.b)

    @staticmethod
    def log(z: 'Dual') -> 'Dual':
        """
        Calculates the logarithm of a dual number.
        
            log(a + b * ϵ) = log(a) + (b / a) * ϵ

        Which is only defined when a > 0.
        
        Parameters:
            z (Dual): The dual number

        Returns:
            Dual: The logarithm of the dual number
        """
        if z.a <= 0:
            raise ValueError("Log undefined: requires a > 0")
        else:
            return Dual(log(z.a), z.b / z.a)
    
    @staticmethod
    def sqrt(z: 'Dual') -> 'Dual':
        """
        Calculates the square root of a dual number.

            sqrt(a + b * ϵ) = sqrt(a) + (b / (2 * sqrt(a))) * ϵ
        
        Which is only defined when a > 0.

        Parameters:
            z (Dual): The dual number
        
        Returns:
            Dual: The square root of the dual number
        """
        if z.a <= 0:
            raise ValueError("Sqrt undefined: requires a > 0")
        else:
            return Dual(sqrt(z.a), z.b / (2 * sqrt(z.a)))
    
    @staticmethod
    def root(n: int, z: 'Dual') -> 'Dual':
        """
        Calculates the n-th root of a dual number.

            root(n, z) = exp((1 / n) * log(z))
        
        Which is only defined when z = 0 or a > 0.

        Expects n > 0.
        
        Parameters:
            z (Dual): The dual number
        
        Returns:
            Dual: The n-th root of the dual number
        """
        if z == 0:
            return Dual(0, 0)
        
        if n <= 0:
            raise ValueError("Root: expected n > 0")
        elif n == 1:
            return z
        else:
            if z.a <= 0:
                raise ValueError("Root undefined: requires a > 0")
            else:
                return Dual.exp((1 / n) * Dual.log(z))

class SplitComplex:
    """
    A class representing a split-complex number in its canonical basis.
    See SplitComplexIdempotent for the idempotent (diagonal) basis.

    Also known as hyperbolic number, perplex number or double number.
    
    A split-complex number (canonical basis) is a number of the form
        a + b * j
    with a, b in R and j^2 = 1.

    Attributes:
        a (float): real component
        b (float): imaginary component
    """
    
    def __init__(self: 'SplitComplex', a: float, b: float) -> 'SplitComplex':
        """
        Initializes a new split-complex number with the given components.

        Parameters:
            a (float): real component
            b (float): imaginary component
        
        Returns:
            SplitComplex: The new split-complex number
        """
        self.a = a
        self.b = b
    
    def __str__(self: 'SplitComplex') -> str:
        return f"({self.a}, {self.b})"
    
    def __repr__(self: 'SplitComplex') -> str:
        return f"SplitComplex({self.a}, {self.b})"
    
    def __eq__(self: 'SplitComplex', rhs) -> bool:
        if isinstance(rhs, (int, float)):
            return (self.a == rhs and self.b == 0)
        elif isinstance(rhs, SplitComplex):
            return (self.a == rhs.a and self.b == rhs.b)
        elif isinstance(rhs, SplitComplexIdempotent):
            zrhs = rhs.to_splitcomplex()
            return (self.a == zrhs.a and self.b == zrhs.b)
        else:
            return NotImplemented

    def __ne__(self: 'SplitComplex', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq
    
    def __abs__(self: 'SplitComplex') -> float:
        return self.norm()
    
    def __neg__(self: 'SplitComplex') -> 'SplitComplex':
        return SplitComplex(-self.a, -self.b)
    
    def __add__(self: 'SplitComplex', rhs) -> 'SplitComplex':
        if isinstance(rhs, (int, float)):
            return SplitComplex(self.a + rhs, self.b)
        elif isinstance(rhs, SplitComplex):
            return SplitComplex(self.a + rhs.a, self.b + rhs.b)
        elif isinstance(rhs, SplitComplexIdempotent):
            zrhs = rhs.to_splitcomplex()
            return SplitComplex(self.a + zrhs.a, self.b + zrhs.b)
        else:
            return NotImplemented
        
    def __radd__(self: 'SplitComplex', lhs) -> 'SplitComplex':
        if isinstance(lhs, (int, float)):
            return SplitComplex(lhs + self.a, self.b)
        elif isinstance(lhs, SplitComplex):
            return SplitComplex(lhs.a + self.a, lhs.b + self.b)
        elif isinstance(lhs, SplitComplexIdempotent):
            zlhs = lhs.to_splitcomplex()
            return SplitComplex(zlhs.a + self.a, zlhs.b + self.b)
        else:
            return NotImplemented
    
    def __sub__(self: 'SplitComplex', rhs) -> 'SplitComplex':
        if isinstance(rhs, (int, float)):
            return SplitComplex(self.a - rhs, self.b)
        elif isinstance(rhs, SplitComplex):
            return SplitComplex(self.a - rhs.a, self.b - rhs.b)
        elif isinstance(rhs, SplitComplexIdempotent):
            zrhs = rhs.to_splitcomplex()
            return SplitComplex(self.a - zrhs.a, self.b - zrhs.b)
        else:
            return NotImplemented
    
    def __rsub__(self: 'SplitComplex', lhs) -> 'SplitComplex':
        if isinstance(lhs, (int, float)):
            return SplitComplex(lhs - self.a, -self.b)
        elif isinstance(lhs, SplitComplex):
            return SplitComplex(lhs.a - self.a, lhs.b - self.b)
        elif isinstance(lhs, SplitComplexIdempotent):
            zlhs = lhs.to_splitcomplex()
            return SplitComplex(zlhs.a - self.a, zlhs.b - self.b)
        else:
            return NotImplemented
    
    def __mul__(self: 'SplitComplex', rhs) -> 'SplitComplex':
        if isinstance(rhs, (int, float)):
            return SplitComplex(self.a * rhs, self.b * rhs)
        elif isinstance(rhs, SplitComplex):
            return SplitComplex(self.a * rhs.a + self.b * rhs.b, self.a * rhs.b + rhs.a * self.b)
        elif isinstance(rhs, SplitComplexIdempotent):
            zrhs = rhs.to_splitcomplex()
            return SplitComplex(self.a * zrhs.a + self.b * zrhs.b, self.a * zrhs.b + zrhs.a * self.b)
        else:
            return NotImplemented

    def __rmul__(self: 'SplitComplex', lhs) -> 'SplitComplex':
        if isinstance(lhs, (int, float)):
            return SplitComplex(lhs * self.a, lhs * self.b)
        elif isinstance(lhs, SplitComplex):
            return SplitComplex(lhs.a * self.a + lhs.b * self.b, lhs.a * self.b + self.a * lhs.b)
        elif isinstance(lhs, SplitComplexIdempotent):
            zlhs = lhs.to_splitcomplex()
            return SplitComplex(zlhs.a * self.a + zlhs.b * self.b, zlhs.a * self.b + self.a * zlhs.b)
        else:
            return NotImplemented
    
    def __truediv__(self: 'SplitComplex', rhs) -> 'SplitComplex':
        if isinstance(rhs, (int, float)):
            return SplitComplex(self.a / rhs, self.b / rhs)
        elif isinstance(rhs, SplitComplex):
            return SplitComplex((self.a * rhs.a - self.b * rhs.b) / (rhs.a**2 - rhs.b**2), (self.b * rhs.a - self.a * rhs.b) / (rhs.a**2 - rhs.b**2))
        elif isinstance(rhs, SplitComplexIdempotent):
            zrhs = rhs.to_splitcomplex()
            return SplitComplex((self.a * zrhs.a - self.b * zrhs.b) / (zrhs.a**2 - zrhs.b**2), (self.b * zrhs.a - self.a * zrhs.b) / (zrhs.a**2 - zrhs.b**2))
        else:
            return NotImplemented
    
    def __rtruediv__(self: 'SplitComplex', lhs) -> 'SplitComplex':
        if isinstance(lhs, (int, float)):
            return SplitComplex((lhs * self.a) / (self.a**2 - self.b**2), (lhs * -self.b) / (self.a**2 - self.b**2))
        elif isinstance(lhs, SplitComplex):
            return SplitComplex((lhs.a * self.a - lhs.b * self.b) / (self.a**2 - self.b**2), (lhs.b * self.a - lhs.a * self.b) / (self.a**2 - self.b**2))
        elif isinstance(lhs, SplitComplexIdempotent):
            zlhs = lhs.to_splitcomplex()
            return SplitComplex((zlhs.a * self.a - zlhs.b * self.b) / (self.a**2 - self.b**2), (zlhs.b * self.a - zlhs.a * self.b) / (self.a**2 - self.b**2))
        else:
            return NotImplemented
        
    def __pow__(self: 'SplitComplex', exp) -> 'SplitComplex':
        if isinstance(exp, int):
            z = SplitComplex(1, 0)

            if exp == 0:
                return z
            
            mul = self if exp > 0 else self.inverse()

            for _ in range(abs(exp)):
                z *= mul
            return z
        elif isinstance(exp, SplitComplex):
            if self.a <= abs(self.b):
                raise ValueError("Pow undefined: base split complex number requires a > |b|")
            else:
                return SplitComplex.exp(exp * SplitComplex.log(self))
        elif isinstance(exp, SplitComplexIdempotent):
            zexp = exp.to_splitcomplex()
            return self**zexp
        else:
            return NotImplemented
    
    def __rpow__(self: 'SplitComplex', base) -> 'SplitComplex':
        if isinstance(base, (int, float)):
            if base <= 0:
                raise ValueError("Pow undefined: base requires x > 0")
            else:
                z = SplitComplex(base, 0)
                return z**self
        elif isinstance(base, SplitComplex):
            return base**self
        else:
            return NotImplemented
    
    def conj(self: 'SplitComplex') -> 'SplitComplex':
        """
        Returns the conjugate of this split-complex number.

        The conjugate of a split-complex number
            a + b * j
        is
            a - b * j

        Returns:
            SplitComplex: The conjugate of this split-complex number
        """
        return SplitComplex(self.a, -self.b)
    
    def norm(self: 'SplitComplex') -> float:
        """
        Calculates the norm of this split-complex number.
        
            a^2 - b^2

        Returns:
            float: The norm of this split-complex number
        """
        return self.a**2 - self.b**2

    def inverse(self: 'SplitComplex') -> 'SplitComplex':
        """
        Calculates the inverse of this split-complex number.

            a / (a^2 - b^2) - b / (a^2 - b^2) * j
        
        Returns:
            SplitComplex: The inverse of this split-complex number
        """
        return SplitComplex(self.a / (self.a**2 - self.b**2), -self.b / (self.a**2 - self.b**2))
    
    def to_idempotent(self: 'SplitComplex') -> 'SplitComplexIdempotent':
        """
        Converts this split-complex number from the canonical basis to the idempotent basis.

            a + b * j
        becomes
            x * e + y * e*
        
        with
            x = a - b
        and
            y = a + b

        Returns:
            SplitComplexIdempotent: The same split-complex number but in the idempotent basis
        """
        return SplitComplexIdempotent(self.a - self.b, self.a + self.b)

    @staticmethod
    def J() -> 'SplitComplex':
        """
        Returns
            0 + 1 * j

        Returns:
            SplitComplex: j
        """
        return SplitComplex(0, 1)

    @staticmethod
    def expJ(x: float) -> 'SplitComplex':
        """
        Calculates e to the power of j * x.

            e^(j * x) = cosh(x) + sinh(x) * j
        
        Parameters:
            x (float): coefficient of j
        
        Returns:
            SplitComplex: The split-complex number representing e^(j * x)
        """
        return SplitComplex(cosh(x), sinh(x))

    @staticmethod
    def exp(z: 'SplitComplex') -> 'SplitComplex':
        """
        Calculates e to the power of a split-complex number.

            e^(a + b * j) = e^a * (cosh(b) + sinh(b) * j)
        
        Parameters:
            z (SplitComplex): The split-complex number in its canonical basis
        
        Returns:
            SplitComplex: e raised to the power of the split-complex number
        """
        return exp(z.a) * SplitComplex(cosh(z.b), sinh(z.b))

    @staticmethod
    def log(z: 'SplitComplex') -> 'SplitComplex':
        """
        Calculates the logarithm of a split-complex number.
        
            log(a + b * j) = log(sqrt(a^2 - b^2)) + atanh(b / a) * j

        Which is only defined when a > |b|.
        
        Parameters:
            z (SplitComplex): The split-complex number in its canonical basis

        Returns:
            SplitComplex: The logarithm of the split-complex number
        """
        if z.a <= abs(z.b):
            raise ValueError("Log undefined: requires a > |b|")
        else:
            return SplitComplex(log(sqrt(z.a**2 - z.b**2)), atanh(z.b / z.a))
    
    @staticmethod
    def sqrt(z: 'SplitComplex') -> 'SplitComplex':
        """
        Calculates the square root of a split-complex number.

            sqrt(a + b * j) = sqrt((a + sqrt(|z|)) / 2) + (b / (2 * sqrt((a + sqrt(|z|)) / 2))) * j

        Which is only defined when a > |b|.
        
        Parameters:
            z (SplitComplex): The split-complex number in its canonical basis
        
        Returns:
            SplitComplex: The square root of the split-complex number
        """
        if z.norm() <= 0:
            raise ValueError("Sqrt undefined: requires |z| > 0")
        else:
            a = sqrt((z.a + sqrt(z.norm())) / 2)
            return SplitComplex(a, z.b / (2 * a))
        
    @staticmethod
    def root(n: int, z: 'SplitComplex') -> 'SplitComplex':
        """
        Calculates the n-th root of a split-complex number.

            root(n, z) = exp((1 / n) * log(z))
        
        Which is only defined when z = 0 or a > |b|.

        Expects n > 0.
        
        Parameters:
            z (SplitComplex): The split-complex number in its canonical basis
        
        Returns:
            SplitComplex: The n-th root of the split-complex number
        """
        if z == 0:
            return SplitComplex(0, 0)
        
        if n <= 0:
            raise ValueError("Root: expected n > 0")
        elif n == 1:
            return z
        else:
            if z.a <= abs(z.b):
                raise ValueError("Root undefined: requires a > |b|")
            else:
                return SplitComplex.exp((1 / n) * SplitComplex.log(z))

    @staticmethod
    def sinh(z: 'SplitComplex') -> 'SplitComplex':
        """
        Calculates sinh of a split-complex number.

            (e^z - e^(-z)) / 2
                
        Parameters:
            z (SplitComplex): The split-complex number in its canonical basis
        
        Returns:
            SplitComplex: sinh of the split-complex number
        """
        return (SplitComplex.exp(z) - SplitComplex.exp(-z)) / 2
    
    @staticmethod
    def cosh(z: 'SplitComplex') -> 'SplitComplex':
        """
        Calculates cosh of a split-complex number.

            (e^z + e^(-z)) / 2
                
        Parameters:
            z (SplitComplex): The split-complex number in its canonical basis
        
        Returns:
            SplitComplex: cosh of the split-complex number
        """
        return (SplitComplex.exp(z) + SplitComplex.exp(-z)) / 2

class SplitComplexIdempotent:
    """
    A class representing a split-complex number in its idempotent (diagonal) basis.
    See SplitComplex for the canonical basis.

    The idempotent basis are:
        e  = 1/2 * (1 - j)
        e* = 1/2 * (1 + j)
        
    These satisfy the properties:
        e^2 = e,    e*^2 = e*,    e * e* = 0,    e + e* = 1

    A split-complex number (idempotent basis) is a number of the form
        a * e + b * e*
    with a, b in R.

    Attributes:
        a (float): coefficient of e
        b (float): coefficient of e*
    """

    def __init__(self: 'SplitComplexIdempotent', a: float, b: float) -> 'SplitComplexIdempotent':
        """
        Initializes a new split-complex number with the given coefficients.

        Parameters:
            a (float): coefficient of e
            b (float): coefficient of e*
        
        Returns:
            SplitComplexIdempotent: The new split-complex number
        """
        self.a = a
        self.b = b
    
    def __str__(self: 'SplitComplexIdempotent') -> str:
        return f"({self.a}, {self.b})"
    
    def __repr__(self: 'SplitComplexIdempotent') -> str:
        return f"SplitComplexIdempotent({self.a}, {self.b})"
    
    def __eq__(self: 'SplitComplexIdempotent', rhs) -> bool:
        if isinstance(rhs, (int, float)):
            zrhs = SplitComplex(rhs, 0).to_idempotent()
            return (self.a == zrhs.a and self.b == zrhs.b)
        elif isinstance(rhs, SplitComplexIdempotent):
            return (self.a == rhs.a and self.b == rhs.b)
        elif isinstance(rhs, SplitComplex):
            zrhs = rhs.to_idempotent()
            return (self.a == zrhs.a and self.b == zrhs.b)
        else:
            return NotImplemented

    def __ne__(self: 'SplitComplexIdempotent', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq
    
    def __abs__(self: 'SplitComplexIdempotent') -> float:
        return self.norm()
    
    def __neg__(self: 'SplitComplexIdempotent') -> 'SplitComplexIdempotent':
        return SplitComplexIdempotent(-self.a, -self.b)
    
    def __add__(self: 'SplitComplexIdempotent', rhs) -> 'SplitComplexIdempotent':
        if isinstance(rhs, (int, float)):
            return SplitComplexIdempotent(self.a + rhs, self.b + rhs)
        elif isinstance(rhs, SplitComplexIdempotent):
            return SplitComplexIdempotent(self.a + rhs.a, self.b + rhs.b)
        elif isinstance(rhs, SplitComplex):
            zrhs = rhs.to_idempotent()
            return SplitComplexIdempotent(self.a + zrhs.a, self.b + zrhs.b)
        else:
            return NotImplemented
    
    def __radd__(self: 'SplitComplexIdempotent', lhs) -> 'SplitComplexIdempotent':
        if isinstance(lhs, (int, float)):
            return SplitComplexIdempotent(lhs + self.a, lhs + self.b)
        elif isinstance(lhs, SplitComplexIdempotent):
            return SplitComplexIdempotent(lhs.a + self.a, lhs.b + self.b)
        elif isinstance(lhs, SplitComplex):
            zlhs = lhs.to_idempotent()
            return SplitComplexIdempotent(zlhs.a + self.a, zlhs.b + self.b)
        else:
            return NotImplemented
    
    def __sub__(self: 'SplitComplexIdempotent', rhs) -> 'SplitComplexIdempotent':
        if isinstance(rhs, (int, float)):
            return SplitComplexIdempotent(self.a - rhs, self.b - rhs)
        elif isinstance(rhs, SplitComplexIdempotent):
            return SplitComplexIdempotent(self.a - rhs.a, self.b - rhs.b)
        elif isinstance(rhs, SplitComplex):
            zrhs = rhs.to_idempotent()
            return SplitComplexIdempotent(self.a - zrhs.a, self.b - zrhs.b)
        else:
            return NotImplemented
    
    def __rsub__(self: 'SplitComplexIdempotent', lhs) -> 'SplitComplexIdempotent':
        if isinstance(lhs, (int, float)):
            return SplitComplexIdempotent(lhs - self.a, lhs - self.b)
        elif isinstance(lhs, SplitComplexIdempotent):
            return SplitComplexIdempotent(lhs.a - self.a, lhs.b - self.b)
        elif isinstance(lhs, SplitComplex):
            zlhs = lhs.to_idempotent()
            return SplitComplexIdempotent(zlhs.a - self.a, zlhs.b - self.b)
        else:
            return NotImplemented
    
    def __mul__(self: 'SplitComplexIdempotent', rhs) -> 'SplitComplexIdempotent':
        if isinstance(rhs, (int, float)):
            return SplitComplexIdempotent(self.a * rhs, self.b * rhs)
        elif isinstance(rhs, SplitComplexIdempotent):
            return SplitComplexIdempotent(self.a * rhs.a, self.b * rhs.b)
        elif isinstance(rhs, SplitComplex):
            zrhs = rhs.to_idempotent()
            return SplitComplexIdempotent(self.a * zrhs.a, self.b * zrhs.b)
        else:
            return NotImplemented
    
    def __rmul__(self: 'SplitComplexIdempotent', lhs) -> 'SplitComplexIdempotent':
        if isinstance(lhs, (int, float)):
            return SplitComplexIdempotent(lhs * self.a, lhs * self.b)
        elif isinstance(lhs, SplitComplexIdempotent):
            return SplitComplexIdempotent(lhs.a * self.a, lhs.b * self.b)
        elif isinstance(lhs, SplitComplex):
            zlhs = lhs.to_idempotent()
            return SplitComplexIdempotent(zlhs.a * self.a, zlhs.b * self.b)
        else:
            return NotImplemented
    
    def __truediv__(self: 'SplitComplexIdempotent', rhs) -> 'SplitComplexIdempotent':
        if isinstance(rhs, (int, float)):
            return SplitComplexIdempotent(self.a / rhs, self.b / rhs)
        elif isinstance(rhs, SplitComplexIdempotent):
            return SplitComplexIdempotent(self.a / rhs.a, self.b / rhs.b)
        elif isinstance(rhs, SplitComplex):
            zrhs = rhs.to_idempotent()
            return SplitComplexIdempotent(self.a / zrhs.a, self.b / zrhs.b)
        else:
            return NotImplemented
    
    def __rtruediv__(self: 'SplitComplexIdempotent', lhs) -> 'SplitComplexIdempotent':
        if isinstance(lhs, (int, float)):
            return SplitComplexIdempotent(lhs / self.a, lhs / self.b)
        elif isinstance(lhs, SplitComplexIdempotent):
            return SplitComplexIdempotent(lhs.a / self.a, lhs.b / self.b)
        elif isinstance(lhs, SplitComplex):
            zlhs = lhs.to_idempotent()
            return SplitComplexIdempotent(zlhs.a / self.a, zlhs.b / self.b)
        else:
            return NotImplemented
    
    def __pow__(self: 'SplitComplexIdempotent', exp) -> 'SplitComplexIdempotent':
        if isinstance(exp, int):
            z = SplitComplexIdempotent(1, 1)

            if exp == 0:
                return z
            
            mul = self if exp > 0 else self.inverse()

            for _ in range(abs(exp)):
                z *= mul
            return z
        elif isinstance(exp, SplitComplexIdempotent):
            if self.a <= 0 or self.b <= 0:
                raise ValueError("Pow undefined: base split-complex number requires a > 0 and b > 0")
            else:
                return SplitComplexIdempotent.exp(exp * SplitComplexIdempotent.log(self))
        elif isinstance(exp, SplitComplex):
            zexp = exp.to_idempotent()
            return self**zexp
        else:
            return NotImplemented
        
    def __rpow__(self: 'SplitComplexIdempotent', base) -> 'SplitComplexIdempotent':
        if isinstance(base, (int, float)):
            if base <= 0:
                raise ValueError("Pow undefined: base requires x > 0")
            else:
                z = SplitComplexIdempotent(base, base)
                return z**self
        elif isinstance(base, SplitComplexIdempotent):
            return base**self
        else:
            return NotImplemented
    
    def conj(self: 'SplitComplexIdempotent') -> 'SplitComplexIdempotent':
        """
        Returns the conjugate of this split-complex number.

        The conjugate of a split-complex number (idempotent basis)
            a * e + b * e*
        is
            b * e + a * e*

        Returns:
            SplitComplexIdempotent: The conjugate of this split-complex number
        """
        return SplitComplexIdempotent(self.b, self.a)
    
    def norm(self: 'SplitComplexIdempotent') -> float:
        """
        Calculates the norm of this split-complex number.
        
            a * b

        Returns:
            float: The norm of this split-complex number
        """
        return self.a * self.b

    def inverse(self: 'SplitComplexIdempotent') -> 'SplitComplexIdempotent':
        """
        Calculates the inverse of this split-complex number.

            1 / a * e + 1 / b * e*
        
        Returns:
            SplitComplexIdempotent: The inverse of this split-complex number
        """
        return SplitComplexIdempotent(1 / self.a, 1 / self.b)
    
    def to_splitcomplex(self: 'SplitComplexIdempotent') -> 'SplitComplex':
        """
        Converts this split-complex number from the idempotent basis to the canonical basis.
        
            a * e + b * e*
        becomes
            x + y * j
        
        with
            x = 1 / 2 * (a + b)
        and
            y = 1 / 2 * (b - a)
        
        Returns:
            SplitComplex: The same split-complex number but in the canonical basis
        """
        return SplitComplex(.5 * (self.a + self.b), .5 * (self.b - self.a))
    
    @staticmethod
    def E() -> 'SplitComplexIdempotent':
        """
        Returns
            1 * e + 0 * e*

        Returns:
            SplitComplexIdempotent: e
        """
        return SplitComplexIdempotent(1, 0)
    
    @staticmethod
    def E2() -> 'SplitComplexIdempotent':
        """
        Returns
            0 * e + 1 * e*

        Returns:
            SplitComplexIdempotent: e*
        """
        return SplitComplexIdempotent(0, 1)
    
    @staticmethod
    def exp(z: 'SplitComplexIdempotent') -> 'SplitComplexIdempotent':
        """
        Calculates e to the power of a split-complex number.

            e^(a * e + b * e*) = e^a * e + e^b * e*
                
        Parameters:
            z (SplitComplexIdempotent): The split-complex number in its idempotent (diagonal) basis
        
        Returns:
            SplitComplexIdempotent: e raised to the power of the split-complex number
        """
        return SplitComplexIdempotent(exp(z.a), exp(z.b))

    @staticmethod
    def log(z: 'SplitComplexIdempotent') -> 'SplitComplexIdempotent':
        """
        Calculates the logarithm of a split-complex number.
        
            log(a * e + b * e*) = log(a) * e + log(b) * e*

        Which is only defined when a > 0 and b > 0.
        
        Parameters:
            z (SplitComplexIdempotent): The split-complex number in its idempotent (diagonal) basis

        Returns:
            SplitComplexIdempotent: The logarithm of the split-complex number
        """
        if z.a <= 0 or z.b <= 0:
            raise ValueError("Log undefined: requires a > 0 and b > 0")
        else:
            return SplitComplexIdempotent(log(z.a), log(z.b))
    
    @staticmethod
    def sqrt(z: 'SplitComplexIdempotent') -> 'SplitComplexIdempotent':
        """
        Calculates the square root of a split-complex number.

            sqrt(a * e + b * e*) = sqrt(a) * e + sqrt(b) * e*

        Which is only defined when a > 0 and b > 0.
        
        Parameters:
            z (SplitComplexIdempotent): The split-complex number in its idempotent (diagonal) basis
        
        Returns:
            SplitComplexIdempotent: The square root of the split-complex number
        """
        if z.a <= 0 or z.b <= 0:
            raise ValueError("Sqrt undefined: requires a > 0 and b > 0")
        else:
            return SplitComplexIdempotent(sqrt(z.a), sqrt(z.b))
        
    @staticmethod
    def root(n: int, z: 'SplitComplexIdempotent') -> 'SplitComplexIdempotent':
        """
        Calculates the n-th root of a split-complex number.

            root(n, z) = exp((1 / n) * log(z))
        
        Which is only defined when z = 0 or a, b > 0.

        Expects n > 0.
        
        Parameters:
            z (SplitComplexIdempotent): The split-complex number in its idempotent (diagonal) basis
        
        Returns:
            SplitComplexIdempotent: The n-th root of the split-complex number
        """
        if z == 0:
            return SplitComplexIdempotent(0, 0)
        
        if n <= 0:
            raise ValueError("Root: expected n > 0")
        elif n == 1:
            return z
        else:
            if z.a <= 0 or z.b <= 0:
                raise ValueError("Root undefined: requires a > 0 and b > 0")
            else:
                return SplitComplexIdempotent.exp((1 / n) * SplitComplexIdempotent.log(z))

    @staticmethod
    def sinh(z: 'SplitComplexIdempotent') -> 'SplitComplexIdempotent':
        """
        Calculates sinh of a split-complex number.

            sinh(a * e + b * e*) = sinh(a) * e + sinh(b) * e*
                
        Parameters:
            z (SplitComplexIdempotent): The split-complex number in its idempotent (diagonal) basis
        
        Returns:
            SplitComplexIdempotent: sinh of the split-complex number
        """
        return SplitComplexIdempotent(sinh(z.a), sinh(z.b))
    
    @staticmethod
    def cosh(z: 'SplitComplexIdempotent') -> 'SplitComplexIdempotent':
        """
        Calculates cosh of a split-complex number.

            cosh(a * e + b * e*) = cosh(a) * e + cosh(b) * e*
                
        Parameters:
            z (SplitComplexIdempotent): The split-complex number in its idempotent (diagonal) basis
        
        Returns:
            SplitComplexIdempotent: cosh of the split-complex number
        """
        return SplitComplexIdempotent(cosh(z.a), cosh(z.b))

class HyperDual:
    """
    A class representing a hyperdual number with real coefficients.

    A hyperdual number is a number of the form
        a + b * ϵ1 + c * ϵ2 + d * (ϵ1 * ϵ2)
    with a, b, c, d in R and ϵ1^2 = ϵ2^2 = (ϵ1 * ϵ2)^2 = 0 and ϵ1 * ϵ2 != 0.

    Attributes:
        a (float): real component
        b (float): 1st dual component
        c (float): 2nd dual component
        d (float): 3rd dual component
    """

    def __init__(self: 'HyperDual', a: float, b: float, c: float, d: float = 0) -> 'HyperDual':
        """
        Initializes a new hyperdual number with the given components.

        Parameters:
            a (float): real component
            b (float): 1st dual component
            c (float): 2nd dual component
            d (float): 3rd dual component or 0 by default
        
        Returns:
            HyperDual: The new hyperdual number
        """
        self.a = a
        self.b = b
        self.c = c
        self.d = d
    
    def __str__(self: 'HyperDual') -> str:
        return f"({self.a}, {self.b}, {self.c}, {self.d})"
    
    def __repr__(self: 'HyperDual') -> str:
        return f"HyperDual({self.a}, {self.b}, {self.c}, {self.d})"
    
    def __eq__(self: 'HyperDual', rhs) -> bool:
        if isinstance(rhs, (int, float)):
            return (self.a == rhs and self.b == 0 and self.c == 0 and self.d == 0)
        elif isinstance(rhs, HyperDual):
            return (self.a == rhs.a and self.b == rhs.b and self.c == rhs.c and self.d == rhs.d)
        else:
            return NotImplemented

    def __ne__(self: 'HyperDual', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq

    def __abs__(self: 'HyperDual') -> float:
        return self.norm()

    def __neg__(self: 'HyperDual') -> 'HyperDual':
        return HyperDual(-self.a, -self.b, -self.c, -self.d)
    
    def __add__(self: 'HyperDual', rhs) -> 'HyperDual':
        if isinstance(rhs, (int, float)):
            return HyperDual(self.a + rhs, self.b, self.c, self.d)
        elif isinstance(rhs, HyperDual):
            return HyperDual(self.a + rhs.a, self.b + rhs.b, self.c + rhs.c, self.d + rhs.d)
        else:
            return NotImplemented
    
    def __radd__(self: 'HyperDual', lhs) -> 'HyperDual':
        if isinstance(lhs, (int, float)):
            return HyperDual(lhs + self.a, self.b, self.c, self.d)
        elif isinstance(lhs, HyperDual):
            return HyperDual(lhs.a + self.a, lhs.b + self.b, lhs.c + self.c, lhs.d + self.d)
        else:
            return NotImplemented
    
    def __sub__(self: 'HyperDual', rhs) -> 'HyperDual':
        if isinstance(rhs, (int, float)):
            return HyperDual(self.a - rhs, self.b, self.c, self.d)
        elif isinstance(rhs, HyperDual):
            return HyperDual(self.a - rhs.a, self.b - rhs.b, self.c - rhs.c, self.d - rhs.d)
        else:
            return NotImplemented
    
    def __rsub__(self: 'HyperDual', lhs) -> 'HyperDual':
        if isinstance(lhs, (int, float)):
            return HyperDual(lhs - self.a, -self.b, -self.c, -self.d)
        elif isinstance(lhs, HyperDual):
            return HyperDual(lhs.a - self.a, lhs.b - self.b, lhs.c - self.c, lhs.d - self.d)
        else:
            return NotImplemented
    
    def __mul__(self: 'HyperDual', rhs) -> 'HyperDual':
        if isinstance(rhs, (int, float)):
            return HyperDual(self.a * rhs, self.b * rhs, self.c * rhs, self.d * rhs)
        elif isinstance(rhs, HyperDual):
            return HyperDual(
                self.a * rhs.a,
                self.a * rhs.b + self.b * rhs.a,
                self.a * rhs.c + self.c * rhs.a,
                self.a * rhs.d + self.b * rhs.c + self.c * rhs.b + self.d * rhs.a
            )
        else:
            return NotImplemented
    
    def __rmul__(self: 'HyperDual', lhs) -> 'HyperDual':
        if isinstance(lhs, (int, float)):
            return HyperDual(lhs * self.a, lhs * self.b, lhs * self.c, lhs * self.d)
        elif isinstance(lhs, HyperDual):
            return HyperDual(
                lhs.a * self.a,
                lhs.a * self.b + lhs.b * self.a,
                lhs.a * self.c + lhs.c * self.a,
                lhs.a * self.d + lhs.b * self.c + lhs.c * self.b + lhs.d * self.a
            )
        else:
            return NotImplemented
        
    def __truediv__(self: 'HyperDual', rhs) -> 'HyperDual':
        if isinstance(rhs, (int, float)):
            return HyperDual(self.a / rhs, self.b / rhs, self.c / rhs, self.d / rhs)
        elif isinstance(rhs, HyperDual):
            return HyperDual(
                self.a / rhs.a,
                (rhs.a * self.b - rhs.b * self.a) / rhs.a**2,
                (rhs.a * self.c - rhs.c * self.a) / rhs.a**2,
                (rhs.a * self.d - rhs.b * self.c - rhs.c * self.b - rhs.d * self.a) / rhs.a**2
            )
        else:
            return NotImplemented
    
    def __rtruediv__(self: 'HyperDual', lhs) -> 'HyperDual':
        if isinstance(lhs, (int, float)):
            return HyperDual(
                lhs / self.a,
                -(lhs * self.b) / self.a**2,
                -(lhs * self.c) / self.a**2,
                (-(lhs * self.d) / self.a**2) + ((2 * lhs * self.b * self.c) / self.a**3)
            )
        elif isinstance(lhs, HyperDual):
            return HyperDual(
                lhs.a / self.a,
                (self.a * lhs.b - self.b * lhs.a) / self.a**2,
                (self.a * lhs.c - self.c * lhs.a) / self.a**2,
                (self.a * lhs.d - self.b * lhs.c - self.c * lhs.b - self.d * lhs.a) / self.a**2
            )
        else:
            return NotImplemented

    def __pow__(self: 'HyperDual', exp) -> 'HyperDual':
        if isinstance(exp, int):
            z = HyperDual(1, 0, 0)

            if exp == 0:
                return z
            
            mul = self if exp > 0 else self.inverse()

            for _ in range(abs(exp)):
                z *= mul
            return z
        elif isinstance(exp, HyperDual):
            if self.a <= 0:
                raise ValueError("Pow undefined: base hyperdual number requires a > 0")
            else:
                return HyperDual.exp(exp * HyperDual.log(self))
        else:
            return NotImplemented

    def __rpow__(self: 'HyperDual', base) -> 'HyperDual':
        if isinstance(base, (int, float)):
            if base <= 0:
                raise ValueError("Pow undefined: base requires x > 0")
            else:
                z = HyperDual(base, 0, 0)
                return z**self
        elif isinstance(base, HyperDual):
            return base**self
        else:
            return NotImplemented

    def conj(self: 'HyperDual') -> 'HyperDual':
        """
        Returns the conjugate of this hyperdual number.

        The conjugate of a hyperdual number
            a + b * ϵ1 + c * ϵ2 + d * (ϵ1 * ϵ2)
        is
            a - b * ϵ1 - c * ϵ2 - d * (ϵ1 * ϵ2)

        Returns:
            HyperDual: The conjugate of this hyperdual number
        """
        return HyperDual(self.a, -self.b, -self.c, -self.d)

    def norm(self: 'HyperDual') -> float:
        """
        Calculates the norm of this hyperdual number.
        
            a^2

        Returns:
            float: The norm of this hyperdual number
        """
        return self.a**2

    def inverse(self: 'HyperDual') -> 'HyperDual':
        """
        Calculates the inverse of this hyperdual number.

            z* / |z|
        
        Returns:
            HyperDual: The inverse of this hyperdual number
        """
        return self.conj() / self.norm()

    @staticmethod
    def E1() -> 'HyperDual':
        """
        Returns
            0 + 1 * ϵ1 + 0 * ϵ2 + 0 * (ϵ1 * ϵ2)

        Returns:
            HyperDual: ϵ1
        """
        return HyperDual(0, 1, 0)
    
    @staticmethod
    def E2() -> 'HyperDual':
        """
        Returns
            0 + 0 * ϵ1 + 1 * ϵ2 + 0 * (ϵ1 * ϵ2)

        Returns:
            HyperDual: ϵ2
        """
        return HyperDual(0, 0, 1)
    
    @staticmethod
    def E1E2() -> 'HyperDual':
        """
        Returns
            0 + 0 * ϵ1 + 0 * ϵ2 + 1 * (ϵ1 * ϵ2)

        Returns:
            HyperDual: (ϵ1 * ϵ2)
        """
        return HyperDual(0, 0, 0, 1)
    
    @staticmethod
    def exp(z: 'HyperDual') -> 'HyperDual':
        """
        Calculates e to the power of a hyperdual number.

            e^(a + b * ϵ1 + c * ϵ2 + d * (ϵ1 * ϵ2)) = e^a * (1 + b * ϵ1 + c * ϵ2 + (d + b * c) * (ϵ1 * ϵ2))
        
        Parameters:
            z (HyperDual): The hyperdual number
        
        Returns:
            HyperDual: e raised to the power of the hyperdual number
        """
        return exp(z.a) * HyperDual(1, z.b, z.c, z.d + z.b * z.c)
    
    @staticmethod
    def log(z: 'HyperDual') -> 'HyperDual':
        """
        Calculates the logarithm of a hyperdual number.
        
            log(a + b * ϵ1 + c * ϵ2 + d * (ϵ1 * ϵ2)) = log(a) + (b / a) * ϵ1 + (c / a) * ϵ2 + (d / a - (b * c) / a^2) * (ϵ1 * ϵ2)

        Which is only defined when a > 0.
        
        Parameters:
            z (HyperDual): The hyperdual number

        Returns:
            HyperDual: The logarithm of the hyperdual number
        """
        if z.a <= 0:
            raise ValueError("Log undefined: requires a > 0")
        else:
            return HyperDual(log(z.a), z.b / z.a, z.c / z.a, z.d / z.a - (z.b * z.c) / z.a**2)
    
    @staticmethod
    def sqrt(z: 'HyperDual') -> 'HyperDual':
        """
        Calculates the square root of a hyperdual number.

            sqrt(a + b * ϵ1 + c * ϵ2 + d * (ϵ1 * ϵ2)) = sqrt(a) + (b / (2 * sqrt(a))) * ϵ1 + (c / (2 * sqrt(a))) * ϵ2 + ((d / (2 * sqrt(a))) - (b * c) / (4 * sqrt(a^3))) * (ϵ1 * ϵ2)
        
        Which is only defined when a > 0.

        Parameters:
            z (HyperDual): The hyperdual number
        
        Returns:
            HyperDual: The square root of the hyperdual number
        """
        if z.a <= 0:
            raise ValueError("Sqrt undefined: requires a > 0")
        else:
            return HyperDual(sqrt(z.a), z.b / (2 * sqrt(z.a)), z.c / (2 * sqrt(z.a)), z.d / (2 * sqrt(z.a)) - (z.b * z.c) / (4 * sqrt(z.a**3)))

    @staticmethod
    def root(n: int, z: 'HyperDual') -> 'HyperDual':
        """
        Calculates the n-th root of a hyperdual number.

            root(n, z) = exp((1 / n) * log(z))
        
        Which is only defined when z = 0 or a > 0.

        Expects n > 0.
        
        Parameters:
            z (HyperDual): The hyperdual number
        
        Returns:
            HyperDual: The n-th root of the hyperdual number
        """
        if z == 0:
            return HyperDual(0, 0, 0)
        
        if n <= 0:
            raise ValueError("Root: expected n > 0")
        elif n == 1:
            return z
        else:
            if z.a <= 0:
                raise ValueError("Root undefined: requires a > 0")
            else:
                return HyperDual.exp((1 / n) * HyperDual.log(z))

class HyperDualComplex:
    """
    A class representing a hyperdual number with complex coefficients.

    A hyperdual number is a number of the form
        a + b * ϵ1 + c * ϵ2 + d * (ϵ1 * ϵ2)
    with a, b, c, d in C and ϵ1^2 = ϵ2^2 = (ϵ1 * ϵ2)^2 = 0 and ϵ1 * ϵ2 != 0.

    Attributes:
        a (complex): real component
        b (complex): 1st dual component
        c (complex): 2nd dual component
        d (complex): 3rd dual component
    """
    
    def __init__(self: 'HyperDualComplex', a: complex, b: complex, c: complex, d: complex = 0) -> 'HyperDualComplex':
        """
        Initializes a new hyperdual number with the given components.

        Parameters:
            a (complex): real component
            b (complex): 1st dual component
            c (complex): 2nd dual component
            d (complex): 3rd dual component or 0 by default
        
        Returns:
            HyperDualComplex: The new hyperdual number
        """
        self.a = a
        self.b = b
        self.c = c
        self.d = d
    
    def __str__(self: 'HyperDualComplex') -> str:
        return f"({self.a}, {self.b}, {self.c}, {self.d})"
    
    def __repr__(self: 'HyperDualComplex') -> str:
        return f"HyperDualComplex({self.a}, {self.b}, {self.c}, {self.d})"

    def __eq__(self: 'HyperDualComplex', rhs) -> bool:
        if isinstance(rhs, (int, float)):
            return (self.a.real == rhs and self.a.imag == 0 and self.b == 0 and self.c == 0 and self.d == 0)
        elif isinstance(rhs, complex):
            return (self.a == rhs and self.b == 0 and self.c == 0 and self.d == 0)
        elif isinstance(rhs, HyperDualComplex):
            return (self.a == rhs.a and self.b == rhs.b and self.c == rhs.c and self.d == rhs.d)
        else:
            return NotImplemented
    
    def __ne__(self: 'HyperDualComplex', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq
    
    def __abs__(self: 'HyperDualComplex') -> float:
        return self.norm()

    def __neg__(self: 'HyperDualComplex') -> 'HyperDualComplex':
        return HyperDualComplex(-self.a, -self.b, -self.c, -self.d)
    
    def __add__(self: 'HyperDualComplex', rhs) -> 'HyperDualComplex':
        if isinstance(rhs, (int, float, complex)):
            return HyperDualComplex(self.a + rhs, self.b, self.c, self.d)
        elif isinstance(rhs, HyperDualComplex):
            return HyperDualComplex(self.a + rhs.a, self.b + rhs.b, self.c + rhs.c, self.d + rhs.d)
        else:
            return NotImplemented
    
    def __radd__(self: 'HyperDualComplex', lhs) -> 'HyperDualComplex':
        if isinstance(lhs, (int, float, complex)):
            return HyperDualComplex(lhs + self.a, self.b, self.c, self.d)
        elif isinstance(lhs, HyperDualComplex):
            return HyperDualComplex(lhs.a + self.a, lhs.b + self.b, lhs.c + self.c, lhs.d + self.d)
        else:
            return NotImplemented
    
    def __sub__(self: 'HyperDualComplex', rhs) -> 'HyperDualComplex':
        if isinstance(rhs, (int, float, complex)):
            return HyperDualComplex(self.a - rhs, self.b, self.c, self.d)
        elif isinstance(rhs, HyperDualComplex):
            return HyperDualComplex(self.a - rhs.a, self.b - rhs.b, self.c - rhs.c, self.d - rhs.d)
        else:
            return NotImplemented
    
    def __rsub__(self: 'HyperDualComplex', lhs) -> 'HyperDualComplex':
        if isinstance(lhs, (int, float, complex)):
            return HyperDualComplex(lhs - self.a, -self.b, -self.c, -self.d)
        elif isinstance(lhs, HyperDualComplex):
            return HyperDualComplex(lhs.a - self.a, lhs.b - self.b, lhs.c - self.c, lhs.d - self.d)
        else:
            return NotImplemented
    
    def __mul__(self: 'HyperDualComplex', rhs) -> 'HyperDualComplex':
        if isinstance(rhs, (int, float, complex)):
            return HyperDualComplex(self.a * rhs, self.b * rhs, self.c * rhs, self.d * rhs)
        elif isinstance(rhs, HyperDualComplex):
            return HyperDualComplex(
                self.a * rhs.a,
                self.a * rhs.b + self.b * rhs.a,
                self.a * rhs.c + self.c * rhs.a,
                self.a * rhs.d + self.b * rhs.c + self.c * rhs.b + self.d * rhs.a
            )
        else:
            return NotImplemented
    
    def __rmul__(self: 'HyperDualComplex', lhs) -> 'HyperDualComplex':
        if isinstance(lhs, (int, float, complex)):
            return HyperDualComplex(lhs * self.a, lhs * self.b, lhs * self.c, lhs * self.d)
        elif isinstance(lhs, HyperDualComplex):
            return HyperDualComplex(
                lhs.a * self.a,
                lhs.a * self.b + lhs.b * self.a,
                lhs.a * self.c + lhs.c * self.a,
                lhs.a * self.d + lhs.b * self.c + lhs.c * self.b + lhs.d * self.a
            )
        else:
            return NotImplemented
    
    def __truediv__(self: 'HyperDualComplex', rhs) -> 'HyperDualComplex':
        if isinstance(rhs, (int, float, complex)):
            return HyperDualComplex(self.a / rhs, self.b / rhs, self.c / rhs, self.d / rhs)
        elif isinstance(rhs, HyperDualComplex):
            return HyperDualComplex(
                self.a / rhs.a,
                (rhs.a * self.b - rhs.b * self.a) / rhs.a**2,
                (rhs.a * self.c - rhs.c * self.a) / rhs.a**2,
                (rhs.a * self.d - rhs.b * self.c - rhs.c * self.b - rhs.d * self.a) / rhs.a**2
            )
        else:
            return NotImplemented
    
    def __rtruediv__(self: 'HyperDualComplex', lhs) -> 'HyperDualComplex':
        if isinstance(lhs, (int, float, complex)):
            return HyperDualComplex(
                lhs / self.a,
                -(lhs * self.b) / self.a**2,
                -(lhs * self.c) / self.a**2,
                (-(lhs * self.d) / self.a**2) + ((2 * lhs * self.b * self.c) / self.a**3)
            )
        elif isinstance(lhs, HyperDualComplex):
            return HyperDualComplex(
                lhs.a / self.a,
                (self.a * lhs.b - self.b * lhs.a) / self.a**2,
                (self.a * lhs.c - self.c * lhs.a) / self.a**2,
                (self.a * lhs.d - self.b * lhs.c - self.c * lhs.b - self.d * lhs.a) / self.a**2
            )
        else:
            return NotImplemented
    
    def __pow__(self: 'HyperDualComplex', exp) -> 'HyperDualComplex':
        if isinstance(exp, int):
            z = HyperDualComplex(1, 0, 0)

            if exp == 0:
                return z
            
            mul = self if exp > 0 else self.inverse()

            for _ in range(abs(exp)):
                z *= mul
            return z
        elif isinstance(exp, HyperDualComplex):
            if self.a == 0:
                raise ValueError("Pow undefined: base hyperdual number requires a != 0")
            else:
                return HyperDualComplex.exp(exp * HyperDualComplex.log(self))
        else:
            return NotImplemented
    
    def __rpow__(self: 'HyperDualComplex', base) -> 'HyperDualComplex':
        if isinstance(base, (int, float, complex)):
            if base == 0:
                raise ValueError("Pow undefined: base requires x != 0")
            else:
                z = HyperDualComplex(base, 0, 0)
                return z**self
        elif isinstance(base, HyperDualComplex):
            return base**self
        else:
            return NotImplemented
    
    def conj(self: 'HyperDualComplex') -> 'HyperDualComplex':
        """
        Returns the conjugate of this hyperdual number.

        The conjugate of a hyperdual number
            a + b * ϵ1 + c * ϵ2 + d * (ϵ1 * ϵ2)
        is
            a - b * ϵ1 - c * ϵ2 - d * (ϵ1 * ϵ2)

        Returns:
            HyperDualComplex: The conjugate of this hyperdual number
        """
        return HyperDualComplex(self.a, -self.b, -self.c, -self.d)
    
    def norm(self: 'HyperDualComplex') -> complex:
        """
        Calculates the norm of this hyperdual number.
        
            a^2

        Returns:
            complex: The norm of this hyperdual number
        """
        return self.a**2
    
    def inverse(self: 'HyperDualComplex') -> 'HyperDualComplex':
        """
        Calculates the inverse of this hyperdual number.

            z* / |z|
        
        Returns:
            HyperDualComplex: The inverse of this hyperdual number
        """
        return self.conj() / self.norm()

    @staticmethod
    def E1() -> 'HyperDualComplex':
        """
        Returns
            0 + 1 * ϵ1 + 0 * ϵ2 + 0 * (ϵ1 * ϵ2)

        Returns:
            HyperDualComplex: ϵ1
        """
        return HyperDualComplex(0, 1, 0)
    
    @staticmethod
    def E2() -> 'HyperDualComplex':
        """
        Returns
            0 + 0 * ϵ1 + 1 * ϵ2 + 0 * (ϵ1 * ϵ2)

        Returns:
            HyperDualComplex: ϵ2
        """
        return HyperDualComplex(0, 0, 1)
    
    @staticmethod
    def E1E2() -> 'HyperDualComplex':
        """
        Returns
            0 + 0 * ϵ1 + 0 * ϵ2 + 1 * (ϵ1 * ϵ2)

        Returns:
            HyperDualComplex: (ϵ1 * ϵ2)
        """
        return HyperDualComplex(0, 0, 0, 1)

    @staticmethod
    def exp(z: 'HyperDualComplex') -> 'HyperDualComplex':
        """
        Calculates e to the power of a hyperdual number.

            e^(a + b * ϵ1 + c * ϵ2 + d * (ϵ1 * ϵ2)) = e^a * (1 + b * ϵ1 + c * ϵ2 + (d + b * c) * (ϵ1 * ϵ2))
        
        Parameters:
            z (HyperDualComplex): The hyperdual number
        
        Returns:
            HyperDualComplex: e raised to the power of the hyperdual number
        """
        return cm.exp(z.a) * HyperDualComplex(1, z.b, z.c, z.d + z.b * z.c)
    
    @staticmethod
    def log(z: 'HyperDualComplex') -> 'HyperDualComplex':
        """
        Calculates the logarithm of a hyperdual number.
        
            log(a + b * ϵ1 + c * ϵ2 + d * (ϵ1 * ϵ2)) = log(a) + (b / a) * ϵ1 + (c / a) * ϵ2 + (d / a - (b * c) / a^2) * (ϵ1 * ϵ2)

        Which is only defined when a != 0.
        
        Parameters:
            z (HyperDualComplex): The hyperdual number

        Returns:
            HyperDualComplex: The logarithm of the hyperdual number
        """
        if z.a == 0:
            raise ValueError("Log undefined: requires a != 0")
        else:
            return HyperDualComplex(cm.log(z.a), z.b / z.a, z.c / z.a, z.d / z.a - (z.b * z.c) / z.a**2)
    
    @staticmethod
    def sqrt(z: 'HyperDualComplex') -> 'HyperDualComplex':
        """
        Calculates the square root of a hyperdual number.

            sqrt(a + b * ϵ1 + c * ϵ2 + d * (ϵ1 * ϵ2)) = sqrt(a) + (b / (2 * sqrt(a))) * ϵ1 + (c / (2 * sqrt(a))) * ϵ2 + ((d / (2 * sqrt(a))) - (b * c) / (4 * sqrt(a^3))) * (ϵ1 * ϵ2)
        
        Which is only defined when a != 0.

        Parameters:
            z (HyperDualComplex): The hyperdual number
        
        Returns:
            HyperDualComplex: The square root of the hyperdual number
        """
        if z.a == 0:
            raise ValueError("Sqrt undefined: requires a != 0")
        else:
            return HyperDual(cm.sqrt(z.a), z.b / (2 * cm.sqrt(z.a)), z.c / (2 * cm.sqrt(z.a)), z.d / (2 * cm.sqrt(z.a)) - (z.b * z.c) / (4 * cm.sqrt(z.a**3)))

    @staticmethod
    def root(n: int, z: 'HyperDualComplex') -> 'HyperDualComplex':
        """
        Calculates the n-th root of a hyperdual number.

            root(n, z) = exp((1 / n) * log(z))
        
        Which is only defined when z = 0 or a != 0.

        Expects n > 0.
        
        Parameters:
            z (HyperDualComplex): The hyperdual number
        
        Returns:
            HyperDualComplex: The n-th root of the hyperdual number
        """
        if z == 0:
            return HyperDualComplex(0, 0, 0)
        
        if n <= 0:
            raise ValueError("Root: expected n > 0")
        elif n == 1:
            return z
        else:
            if z.a == 0:
                raise ValueError("Root undefined: requires a != 0")
            else:
                return HyperDualComplex.exp((1 / n) * HyperDualComplex.log(z))

class Quat:
    """
    A class representing a quaternion.
    
    A quaternion is a number of the form
        w + x * i + y * j + z * k
    with w, x, y, z in R and i^2 = j^2 = k^2 = ijk = -1.

    Attributes:
        w (float): real component
        x (float): 1st imaginary component
        y (float): 2nd imaginary component
        z (float): 3rd imaginary component
    """
    
    def __init__(self: 'Quat', w: float, x: float, y: float, z: float) -> 'Quat':
        """
        Initializes a new quaternion with the given components.

        Parameters:
            w (float): real component
            x (float): 1st imaginary component
            y (float): 2nd imaginary component
            z (float): 3rd imaginary component
        
        Returns:
            Quat: The new quaternion
        """
        self.w = w
        self.x = x
        self.y = y
        self.z = z

    def __str__(self: 'Quat') -> str:
        return f"({self.w}, {self.x}, {self.y}, {self.z})"

    def __repr__(self: 'Quat') -> str:
        return f"Quat({self.w}, {self.x}, {self.y}, {self.z})"
    
    def __eq__(self: 'Quat', rhs) -> bool:
        if isinstance(rhs, (int, float)):
            return (self.w == rhs and self.x == 0 and self.y == 0 and self.z == 0)
        elif isinstance(rhs, Quat):
            return (self.w == rhs.w and self.x == rhs.x and self.y == rhs.y and self.z == rhs.z)
        else:
            return NotImplemented

    def __ne__(self: 'Quat', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq

    def __abs__(self: 'Quat') -> float:
        return self.norm()

    def __neg__(self: 'Quat') -> 'Quat':
        return Quat(-self.w, -self.x, -self.y, -self.z)

    def __add__(self: 'Quat', rhs) -> 'Quat':
        if isinstance(rhs, (int, float)):
            return Quat(self.w + rhs, self.x, self.y, self.z)
        elif isinstance(rhs, Quat):
            return Quat(self.w + rhs.w, self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
        else:
            return NotImplemented
    
    def __radd__(self: 'Quat', lhs) -> 'Quat':
        if isinstance(lhs, (int, float)):
            return Quat(lhs + self.w, self.x, self.y, self.z)
        elif isinstance(lhs, Quat):
            return Quat(lhs.w + self.w, lhs.x + self.x, lhs.y + self.y, lhs.z + self.z)
        else:
            return NotImplemented
    
    def __sub__(self: 'Quat', rhs) -> 'Quat':
        if isinstance(rhs, (int, float)):
            return Quat(self.w - rhs, self.x, self.y, self.z)
        elif isinstance(rhs, Quat):
            return Quat(self.w - rhs.w, self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
        else:
            return NotImplemented

    def __rsub__(self: 'Quat', lhs) -> 'Quat':
        if isinstance(lhs, (int, float)):
            return Quat(lhs - self.w, -self.x, -self.y, -self.z)
        elif isinstance(lhs, Quat):
            return Quat(lhs.w - self.w, lhs.x - self.x, lhs.y - self.y, lhs.z - self.z)
        else:
            return NotImplemented

    def __mul__(self: 'Quat', rhs) -> 'Quat':
        if isinstance(rhs, (int, float)):
            return Quat(self.w * rhs, self.x * rhs, self.y * rhs, self.z * rhs)
        elif isinstance(rhs, Quat):
            return Quat(
                self.w * rhs.w - self.x * rhs.x - self.y * rhs.y - self.z * rhs.z,
                self.w * rhs.x + self.x * rhs.w + self.y * rhs.z - self.z * rhs.y,
                self.w * rhs.y - self.x * rhs.z + self.y * rhs.w + self.z * rhs.x,
                self.w * rhs.z + self.x * rhs.y - self.y * rhs.x + self.z * rhs.w
            )
        else:
            return NotImplemented
    
    def __rmul__(self: 'Quat', lhs) -> 'Quat':
        if isinstance(lhs, (int, float)):
            return Quat(lhs * self.w, lhs * self.x, lhs * self.y, lhs * self.z)
        elif isinstance(lhs, Quat):
            return Quat(
                lhs.w * self.w - lhs.x * self.x - lhs.y * self.y - lhs.z * self.z,
                lhs.w * self.x + lhs.x * self.w + lhs.y * self.z - lhs.z * self.y,
                lhs.w * self.y - lhs.x * self.z + lhs.y * self.w + lhs.z * self.x,
                lhs.w * self.z + lhs.x * self.y - lhs.y * self.x + lhs.z * self.w
            )
        else:
            return NotImplemented
        
    # NOTE: __truediv__ and __rtruediv__ are emitted to force multiplication with the inverse which should aid with clarity

    def __pow__(self: 'Quat', exp) -> 'Quat':
        if isinstance(exp, int):
            q = Quat(1, 0, 0, 0)

            if exp == 0:
                return q
            
            mul = self if exp > 0 else self.inverse()

            for _ in range(abs(exp)):
                q *= mul
            return q
        elif isinstance(exp, Quat):
            if self == 0:
                raise ValueError("Pow undefined: base quaternion requires q != 0")
            else:
                return Quat.exp(exp * Quat.log(self))
        else:
            return NotImplemented
        
    def __rpow__(self: 'Quat', base) -> 'Quat':
        if isinstance(base, (int, float)):
            if base == 0:
                raise ValueError("Pow undefined: base requires x != 0")
            else:
                q = Quat(base, 0, 0, 0)
                return q**self
        elif isinstance(base, Quat):
            return base**self
        else:
            return NotImplemented

    def imag(self: 'Quat') -> Vec3:
        """
        Returns the imaginary vector part as Vec3.

        Returns:
            Vec3: The imaginary vector part
        """
        return Vec3(self.x, self.y, self.z)
    
    def conj(self: 'Quat') -> 'Quat':
        """
        Returns the conjugate of this quaternion.

        The conjugate of a quaternion
            w + x * i + y * j + z * k
        is
            w - x * i - y * j - z * k

        Returns:
            Quat: The conjugate of this quaternion
        """
        return Quat(self.w, -self.x, -self.y, -self.z)
    
    def norm(self: 'Quat') -> float:
        """
        Calculates the norm of this quaternion.
        
            sqrt(w^2 + x^2 + y^2 + z^2)

        Returns:
            float: The norm of this quaternion
        """
        return sqrt(self.w**2 + self.x**2 + self.y**2 + self.z**2)
    
    def norm_imm(self: 'Quat') -> float:
        """
        Calculates the norm of this quaternion's imaginary vector part.
        
            sqrt(x^2 + y^2 + z^2)

        Returns:
            float: The norm of this quaternion's imaginary vector part
        """
        return sqrt(self.x**2 + self.y**2 + self.z**2)

    def inverse(self: 'Quat') -> 'Quat':
        """
        Calculates the inverse of this quaternion.

            q* / |q|^2
        
        Returns:
            Quat: The inverse of this quaternion
        """
        return self.conj() * (1 / self.norm()**2)

    @staticmethod
    def I() -> 'Quat':
        """
        Returns
            0 + 1 * i + 0 * j + 0 * k

        Returns:
            Quat: i
        """
        return Quat(0, 1, 0, 0)
    
    @staticmethod
    def J() -> 'Quat':
        """
        Returns
            0 + 0 * i + 1 * j + 0 * k

        Returns:
            Quat: j
        """
        return Quat(0, 0, 1, 0)
    
    @staticmethod
    def K() -> 'Quat':
        """
        Returns
            0 + 0 * i + 0 * j + 1 * k

        Returns:
            Quat: k
        """
        return Quat(0, 0, 0, 1)

    @staticmethod
    def from_vec3(vec: Vec3) -> 'Quat':
        """
        Initializes a quaternion from a Vec3.

        A Vec3
            (v_x, v_y, v_z)
        results in
            v_x * i + v_y * j + v_z * k
        
        Parameters:
            vec (Vec3): The 3-dimensional vector
        
        Returns:
            Quat: An imaginary quaternion with vector part equal to the Vec3
        """
        return Quat(0, vec.x, vec.y, vec.z)
    
    @staticmethod
    def exp(q: 'Quat') -> 'Quat':
        """
        Calculates e to the power of a quaternion (v is the imaginary vector part).

            e^q = e^w * (cos(|v|) + v/|v| * sin(|v|))
        
        Parameters:
            q (Quat): The quaternion
        
        Returns:
            Quat: e raised to the power of the quaternion
        """
        if q.norm_imm() == 0:
            return Quat(exp(q.w), 0, 0, 0)
        else:
            return exp(q.w) * (cos(q.norm_imm()) + Quat(0, q.x, q.y, q.z) * (1 / q.norm_imm()) * sin(q.norm_imm()))

    @staticmethod
    def log(q: 'Quat') -> 'Quat':
        """
        Calculates the logarithm of a quaternion (v is the imaginary vector part).

            log(q) = log(|q|) + v/|v| * atan2(|v|, w)

        Which is only defined when q != 0.
        
        Parameters:
            q (Quat): The quaternion

        Returns:
            Quat: The logarithm of the quaternion
        """
        imm_norm = q.norm_imm()

        if imm_norm == 0:
            if q.w > 0:
                return Quat(log(q.w), 0, 0, 0)
            elif q.w < 0:
                # NOTE: many possible values, choose principal branch with i-direction
                return Quat(log(abs(q.w)), pi, 0, 0)
            else:
                raise ValueError("Log undefined: requires q != 0")
        
        theta = atan2(imm_norm, q.w) # (-pi, pi] principal value
        return Quat(log(q.norm()), q.x * (theta / imm_norm), q.y * (theta / imm_norm), q.z * (theta / imm_norm))

    @staticmethod
    def sqrt(q: 'Quat') -> 'Quat':
        """
        Calculates the square root of a quaternion (v is the imaginary vector part).

            sqrt(q) = sqrt((|q| + w) / 2) + v/|v| * sqrt((|q| - w) / 2)
        
        Parameters:
            q (Quat): The quaternion
        
        Returns:
            Quat: The square root of the quaternion
        """
        imm_norm = q.norm_imm()

        if imm_norm == 0:
            if q.w >= 0:
                return Quat(sqrt(q.w), 0, 0, 0)
            elif q.w < 0:
                # NOTE: many possible values, choose principal branch with i-direction
                return Quat(0, sqrt(abs(q.w)), 0, 0)

        q_norm = q.norm()
        f = sqrt((q_norm - q.w) / 2)
        v = Quat(0, q.x / imm_norm, q.y / imm_norm, q.z / imm_norm)
        return Quat(sqrt((q_norm + q.w) / 2), 0, 0, 0) + f * v

    @staticmethod
    def root(n: int, q: 'Quat') -> 'Quat':
        """
        Calculates the n-th root of a quaternion.

            root(n, q) = exp((1 / n) * log(q))

        Expects n > 0.
        
        Parameters:
            q (Quat): The quaternion
        
        Returns:
            Quat: The n-th root of the quaternion
        """
        if q == 0:
            return Quat(0, 0, 0, 0)

        if n <= 0:
            raise ValueError("Root: expected n > 0")
        elif n == 1:
            return q
        else:
            return Quat.exp((1 / n) * Quat.log(q))

    @staticmethod
    def exp_rot(vec: Vec3, phi: float) -> 'Quat':
        """
        Helper function that creates a rotation quaternion.
        
            e^(v * phi) = cos(phi) + v * sin(phi)

        Parameters:
            vec (Vec3): The rotation axis with length 1 (normalized)
            phi (float): The rotation angle
        
        Returns:
            Quat: The rotation quaternion that rotates around the given axis with the specified angle
        """
        q = sin(phi) * Quat.from_vec3(vec)
        return cos(phi) + q

class SplitQuat:
    """
    A class representing a split-quaternion.
    
    A split-quaternion is a number of the form
        w + x * i + y * j + z * k
    with w, x, y, z in R and i^2 = j^2 = 1 and k^2 = ijk = -1.

    Attributes:
        w (float): real component
        x (float): 1st imaginary component
        y (float): 2nd imaginary component
        z (float): 3rd imaginary component
    """
    
    def __init__(self: 'SplitQuat', w: float, x: float, y: float, z: float) -> 'SplitQuat':
        """
        Initializes a new split-quaternion with the given components.

        Parameters:
            w (float): real component
            x (float): 1st imaginary component
            y (float): 2nd imaginary component
            z (float): 3rd imaginary component
        
        Returns:
            SplitQuat: The new split-quaternion
        """
        self.w = w
        self.x = x
        self.y = y
        self.z = z
    
    def __str__(self: 'SplitQuat') -> str:
        return f"({self.w}, {self.x}, {self.y}, {self.z})"
    
    def __repr__(self: 'SplitQuat') -> str:
        return f"SplitQuat({self.w}, {self.x}, {self.y}, {self.z})"
    
    def __eq__(self: 'SplitQuat', rhs) -> bool:
        if isinstance(rhs, (int, float)):
            return (self.w == rhs and self.x == 0 and self.y == 0 and self.z == 0)
        elif isinstance(rhs, SplitQuat):
            return (self.w == rhs.w and self.x == rhs.x and self.y == rhs.y and self.z == rhs.z)
        else:
            return NotImplemented

    def __ne__(self: 'SplitQuat', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq

    def __abs__(self: 'SplitQuat') -> float:
        return self.norm()
    
    def __neg__(self: 'SplitQuat') -> 'SplitQuat':
        return SplitQuat(-self.w, -self.x, -self.y, -self.z)
    
    def __add__(self: 'SplitQuat', rhs) -> 'SplitQuat':
        if isinstance(rhs, (int, float)):
            return SplitQuat(self.w + rhs, self.x, self.y, self.z)
        elif isinstance(rhs, SplitQuat):
            return SplitQuat(self.w + rhs.w, self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
        else:
            return NotImplemented
    
    def __radd__(self: 'SplitQuat', lhs) -> 'SplitQuat':
        if isinstance(lhs, (int, float)):
            return SplitQuat(lhs + self.w, self.x, self.y, self.z)
        elif isinstance(lhs, SplitQuat):
            return SplitQuat(lhs.w + self.w, lhs.x + self.x, lhs.y + self.y, lhs.z + self.z)
        else:
            return NotImplemented
    
    def __sub__(self: 'SplitQuat', rhs) -> 'SplitQuat':
        if isinstance(rhs, (int, float)):
            return SplitQuat(self.w - rhs, self.x, self.y, self.z)
        elif isinstance(rhs, SplitQuat):
            return SplitQuat(self.w - rhs.w, self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
        else:
            return NotImplemented
    
    def __rsub__(self: 'SplitQuat', lhs) -> 'SplitQuat':
        if isinstance(lhs, (int, float)):
            return SplitQuat(lhs - self.w, -self.x, -self.y, -self.z)
        elif isinstance(lhs, SplitQuat):
            return SplitQuat(lhs.w - self.w, lhs.x - self.x, lhs.y - self.y, lhs.z - self.z)
        else:
            return NotImplemented
    
    def __mul__(self: 'SplitQuat', rhs) -> 'SplitQuat':
        if isinstance(rhs, (int, float)):
            return SplitQuat(self.w * rhs, self.x * rhs, self.y * rhs, self.z * rhs)
        elif isinstance(rhs, SplitQuat):
            return SplitQuat(
                self.w * rhs.w + self.x * rhs.x + self.y * rhs.y - self.z * rhs.z,
                self.w * rhs.x + self.x * rhs.w + self.z * rhs.y - self.y * rhs.z,
                self.w * rhs.y + self.y * rhs.w + self.x * rhs.z - self.z * rhs.x,
                self.w * rhs.z + self.z * rhs.w + self.x * rhs.y - self.y * rhs.x
            )
        else:
            return NotImplemented
    
    def __rmul__(self: 'SplitQuat', lhs) -> 'SplitQuat':
        if isinstance(lhs, (int, float)):
            return SplitQuat(lhs * self.w, lhs * self.x, lhs * self.y, lhs * self.z)
        elif isinstance(lhs, SplitQuat):
            return SplitQuat(
                lhs.w * self.w + lhs.x * self.x + lhs.y * self.y - lhs.z * self.z,
                lhs.w * self.x + lhs.x * self.w + lhs.z * self.y - lhs.y * self.z,
                lhs.w * self.y + lhs.y * self.w + lhs.x * self.z - lhs.z * self.x,
                lhs.w * self.z + lhs.z * self.w + lhs.x * self.y - lhs.y * self.x
            )
        else:
            return NotImplemented

    # NOTE: __truediv__ and __rtruediv__ are emitted to force multiplication with the inverse which should aid with clarity

    def __pow__(self: 'SplitQuat', exp) -> 'SplitQuat':
        if isinstance(exp, int):
            q = SplitQuat(1, 0, 0, 0)

            if exp == 0:
                return q
            
            mul = self if exp > 0 else self.inverse()

            for _ in range(abs(exp)):
                q *= mul
            return q
        elif isinstance(exp, SplitQuat):
            if self == 0 or self.norm() <= 0:
                raise ValueError("Pow undefined: base split-quaternion requires q != 0 and |q| > 0")
            else:
                return SplitQuat.exp(exp * SplitQuat.log(self))
        else:
            return NotImplemented
    
    def __rpow__(self: 'SplitQuat', base) -> 'SplitQuat':
        if isinstance(base, (int, float)):
            if base <= 0:
                raise ValueError("Pow undefined: base requires x > 0")
            else:
                q = SplitQuat(base, 0, 0, 0)
                return q**self
        elif isinstance(base, SplitQuat):
            return base**self
        else:
            return NotImplemented

    def conj(self: 'SplitQuat') -> 'SplitQuat':
        """
        Returns the conjugate of this split-quaternion.

        The conjugate of a split-quaternion
            w + x * i + y * j + z * k
        is
            w - x * i - y * j - z * k

        Returns:
            SplitQuat: The conjugate of this split-quaternion
        """
        return SplitQuat(self.w, -self.x, -self.y, -self.z)

    def norm(self: 'SplitQuat') -> float:
        """
        Calculates the norm of this split-quaternion.
        
            w^2 - x^2 - y^2 + z^2

        Returns:
            float: The norm of this split-quaternion
        """
        return self.w**2 - self.x**2 - self.y**2 + self.z**2
    
    def inverse(self: 'SplitQuat') -> 'SplitQuat':
        """
        Calculates the inverse of this split-quaternion.

            q* / |q|
        
        Returns:
            SplitQuat: The inverse of this split-quaternion
        """
        return self.conj() * (1 / self.norm())

    @staticmethod
    def I() -> 'SplitQuat':
        """
        Returns
            0 + 1 * i + 0 * j + 0 * k

        Returns:
            SplitQuat: i
        """
        return SplitQuat(0, 1, 0, 0)
    
    @staticmethod
    def J() -> 'SplitQuat':
        """
        Returns
            0 + 0 * i + 1 * j + 0 * k

        Returns:
            SplitQuat: j
        """
        return SplitQuat(0, 0, 1, 0)
    
    @staticmethod
    def K() -> 'SplitQuat':
        """
        Returns
            0 + 0 * i + 0 * j + 1 * k

        Returns:
            SplitQuat: k
        """
        return SplitQuat(0, 0, 0, 1)
    
    @staticmethod
    def expI(x: float) -> 'SplitQuat':
        """
        Calculates the exponential in the i-direction.

            e^(x * i) = cosh(x) + sinh(x) * i

        Parameters:
            x (float): coefficient of i

        Returns:
            SplitQuat: The split-quaternion representing e^(x * i)
        """
        return SplitQuat(cosh(x), sinh(x), 0, 0)
    
    @staticmethod
    def expJ(x: float) -> 'SplitQuat':
        """
        Calculates the exponential in the j-direction.

            e^(x * j) = cosh(x) + sinh(x) * j

        Parameters:
            x (float): coefficient of j

        Returns:
            SplitQuat: The split-quaternion representing e^(x * j)
        """
        return SplitQuat(cosh(x), 0, sinh(x), 0)
    
    @staticmethod
    def expK(x: float) -> 'SplitQuat':
        """
        Calculates the exponential in the k-direction.

            e^(x * k) = cos(x) + sin(x) * k

        Parameters:
            x (float): coefficient of k

        Returns:
            SplitQuat: The split-quaternion representing e^(x * k)
        """
        return SplitQuat(cos(x), 0, 0, sin(x))
    
    @staticmethod
    def exp(q: 'SplitQuat') -> 'SplitQuat':
        """
        Calculates e to the power of a split-quaternion (r^2 = x^2 + y^2 - z^2, v is the imaginary vector part).

            e^q = 
                e^w * (cosh(r) + v/r * sinh(r)) when r^2 > 0
                e^w * (cos(r)  + v/r * sin(r)) when r^2 < 0
                e^w * (1 + v) when r^2 = 0

        Parameters:
            q (SplitQuat): The split-quaternion
        
        Returns:
            SplitQuat: e raised to the power of the split-quaternion
        """
        r2 = q.x**2 + q.y**2 - q.z**2
        v = SplitQuat(0, q.x, q.y, q.z)
        if r2 > 0:
            r = sqrt(r2)
            return exp(q.w) * (cosh(r) + v * sinh(r) * (1/r))
        elif r2 < 0:
            r = sqrt(-r2)
            return exp(q.w) * (cos(r) + v * sin(r) * (1/r))
        elif r2 == 0:
            return exp(q.w) * (1 + v)
    
    @staticmethod
    def log(q: 'SplitQuat') -> 'SplitQuat':
        """
        Calculates the logarithm of a split-quaternion (r^2 = x^2 + y^2 - z^2, v is the imaginary vector part).

            log(q) = .5 * log(|q|) +
                v * (atanh(r / w) / r) when r^2 > 0
                v * (atan2(r, w) / r) when r^2 < 0
                v * (1 / w) when r^2 = 0

        Which is only defined when q != 0 and |q| > 0.
        
        Parameters:
            q (SplitQuat): The split-quaternion

        Returns:
            SplitQuat: The logarithm of the split-quaternion
        """
        n = q.norm()
        if n <= 0:
            raise ValueError("Log undefined: requires |q| > 0")
        
        w = .5 * log(n)

        r2 = q.x**2 + q.y**2 - q.z**2
        if r2 > 0:
            r = sqrt(r2)
            phi = atanh(r / q.w)
            return w + SplitQuat(0, q.x, q.y, q.z) * (phi / r)
        elif r2 < 0:
            # NOTE: principal branch as complex log is multivalued
            r = sqrt(-r2)
            theta = atan2(r, q.w)
            return w + SplitQuat(0, q.x, q.y, q.z) * (theta / r)
        elif r2 == 0:
            if q.x == q.y == q.z == 0:
                if q.w > 0:
                    return SplitQuat(log(q.w), 0, 0, 0)
                elif q.w < 0:
                    # NOTE: many possible values, choose principal branch with k-direction
                    return SplitQuat(log(abs(q.w)), 0, 0, pi)
                else:
                    raise ValueError("Log undefined: requires q != 0")
            else:
                return w + SplitQuat(0, q.x, q.y, q.z) * (1 / q.w)
    
    @staticmethod
    def sqrt(q: 'SplitQuat') -> 'SplitQuat':
        """
        Calculates the square root of a split-quaternion.

            sqrt(q) = exp((1 / 2) * log(q))
        
        Which is only defined when q = 0 or |q| > 0.
        
        Parameters:
            q (SplitQuat): The split-quaternion
        
        Returns:
            SplitQuat: The square root of the split-quaternion
        """
        if q == 0:
            return SplitQuat(0, 0, 0, 0)
        
        if q.norm() <= 0:
            raise ValueError("Sqrt undefined: requires |q| > 0")
        else:
            return SplitQuat.exp(.5 * SplitQuat.log(q))
    
    @staticmethod
    def root(n: int, q: 'SplitQuat') -> 'SplitQuat':
        """
        Calculates the n-th root of a split-quaternion.

            root(n, q) = exp((1 / n) * log(q))
        
        Which is only defined when q = 0 or |q| > 0.

        Expects n > 0.
        
        Parameters:
            q (SplitQuat): The split-quaternion
        
        Returns:
            SplitQuat: The n-th root of the split-quaternion
        """
        if q == 0:
            return SplitQuat(0, 0, 0, 0)

        if n <= 0:
            raise ValueError("Root: expected n > 0")
        elif n == 1:
            return q
        else:
            if q.norm() <= 0:
                raise ValueError("Root undefined: requires |q| > 0")
            else:
                return SplitQuat.exp((1 / n) * SplitQuat.log(q))

class DualComplex:
    """
    A class representing a dual complex number.
    
    A dual complex number is a number of the form
        a + b * ϵ
    with a, b in C and ϵ^2 = 0.

    Attributes:
        a (complex): real component
        b (complex): dual component
    """
    
    def __init__(self: 'DualComplex', a: complex, b: complex) -> 'DualComplex':
        """
        Initializes a new dual complex number with the given components.

        Parameters:
            a (complex): real component
            b (complex): dual component
        
        Returns:
            DualComplex: The new dual complex number
        """
        self.a = a
        self.b = b
    
    def __str__(self: 'DualComplex') -> str:
        return f"({self.a}, {self.b})"
    
    def __repr__(self: 'DualComplex') -> str:
        return f"DualComplex({self.a}, {self.b})"
    
    def __eq__(self: 'DualComplex', rhs) -> bool:
        if isinstance(rhs, (int, float)):
            return (self.a.real == rhs and self.a.imag == 0 and self.b == 0)
        elif isinstance(rhs, complex):
            return (self.a == rhs and self.b == 0)
        elif isinstance(rhs, DualComplex):
            return (self.a == rhs.a and self.b == rhs.b)
        else:
            return NotImplemented

    def __ne__(self: 'DualComplex', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq

    def __abs__(self: 'DualComplex') -> 'DualComplex':
        return self.norm()

    def __neg__(self: 'DualComplex') -> 'DualComplex':
        return DualComplex(-self.a, -self.b)
    
    def __add__(self: 'DualComplex', rhs) -> 'DualComplex':
        if isinstance(rhs, (int, float, complex)):
            return DualComplex(self.a + rhs, self.b)
        elif isinstance(rhs, DualComplex):
            return DualComplex(self.a + rhs.a, self.b + rhs.b)
        else:
            return NotImplemented
    
    def __radd__(self: 'DualComplex', lhs) -> 'DualComplex':
        if isinstance(lhs, (int, float, complex)):
            return DualComplex(lhs + self.a, self.b)
        elif isinstance(lhs, DualComplex):
            return DualComplex(lhs.a + self.a, lhs.b + self.b)
        else:
            return NotImplemented
    
    def __sub__(self: 'DualComplex', rhs) -> 'DualComplex':
        if isinstance(rhs, (int, float, complex)):
            return DualComplex(self.a - rhs, self.b)
        elif isinstance(rhs, DualComplex):
            return DualComplex(self.a - rhs.a, self.b - rhs.b)
        else:
            return NotImplemented
    
    def __rsub__(self: 'DualComplex', lhs) -> 'DualComplex':
        if isinstance(lhs, (int, float, complex)):
            return DualComplex(lhs - self.a, -self.b)
        elif isinstance(lhs, DualComplex):
            return DualComplex(lhs.a - self.a, lhs.b - self.b)
        else:
            return NotImplemented
    
    def __mul__(self: 'DualComplex', rhs) -> 'DualComplex':
        if isinstance(rhs, (int, float, complex)):
            return DualComplex(self.a * rhs, self.b * rhs)
        elif isinstance(rhs, DualComplex):
            return DualComplex(self.a * rhs.a, self.a * rhs.b + self.b * rhs.a)
        else:
            return NotImplemented
    
    def __rmul__(self: 'DualComplex', lhs) -> 'DualComplex':
        if isinstance(lhs, (int, float, complex)):
            return DualComplex(lhs * self.a, lhs * self.b)
        elif isinstance(lhs, DualComplex):
            return DualComplex(lhs.a * self.a, lhs.a * self.b + lhs.b * self.a)
        else:
            return NotImplemented
    
    def __truediv__(self: 'DualComplex', rhs) -> 'DualComplex':
        if isinstance(rhs, (int, float, complex)):
            return DualComplex(self.a / rhs, self.b / rhs)
        elif isinstance(rhs, DualComplex):
            return DualComplex(self.a / rhs.a, (self.b * rhs.a - self.a * rhs.b) / rhs.a**2)
        else:
            return NotImplemented
    
    def __rtruediv__(self: 'DualComplex', lhs) -> 'DualComplex':
        if isinstance(lhs, (int, float, complex)):
            return DualComplex(lhs / self.a, (-lhs * self.b) / self.a**2)
        elif isinstance(lhs, DualComplex):
            return DualComplex(lhs.a / self.a, (lhs.b * self.a - lhs.a * self.b) / self.a**2)
        else:
            return NotImplemented

    def __pow__(self: 'DualComplex', exp) -> 'DualComplex':
        if isinstance(exp, int):
            z = DualComplex(1, 0)

            if exp == 0:
                return z
            
            mul = self if exp > 0 else self.inverse()

            for _ in range(abs(exp)):
                z *= mul
            return z
        elif isinstance(exp, DualComplex):
            if self.a == 0:
                raise ValueError("Pow undefined: base dual complex number requires a != 0")
            else:
                return DualComplex.exp(exp * DualComplex.log(self))
        else:
            return NotImplemented
    
    def __rpow__(self: 'DualComplex', base) -> 'DualComplex':
        if isinstance(base, (int, float, complex)):
            if base == 0:
                raise ValueError("Pow undefined: base requires x != 0")
            else:
                z = DualComplex(base, 0)
                return z**self
        elif isinstance(base, DualComplex):
            return base**self
        else:
            return NotImplemented

    def conj_complex(self: 'DualComplex') -> 'DualComplex':
        """
        Returns the complex conjugate of this dual complex number.

        The complex conjugate of a dual complex number
            a + b * ϵ
        is
            a* + b* * ϵ

        Returns:
            DualComplex: The complex conjugate of this dual complex number
        """
        return DualComplex(self.a.conjugate(), self.b.conjugate())
    
    def conj_dual(self: 'DualComplex') -> 'DualComplex':
        """
        Returns the dual conjugate of this dual complex number.

        The dual conjugate of a dual complex number
            a + b * ϵ
        is
            a - b * ϵ

        Returns:
            DualComplex: The dual conjugate of this dual complex number
        """
        return DualComplex(self.a, -self.b)
    
    def conj_both(self: 'DualComplex') -> 'DualComplex':
        """
        Returns the full (complex and dual) conjugate of this dual complex number.

        The full conjugate of a dual complex number
            a + b * ϵ
        is
            a* - b* * ϵ

        Returns:
            DualComplex: The full conjugate of this dual complex number
        """
        return DualComplex(self.a.conjugate(), -self.b.conjugate())

    def norm(self: 'DualComplex') -> 'DualComplex':
        """
        Calculates the norm of this dual complex number.
        
            (a * a*) + (a * b* + b * a*) * ϵ

        Returns:
            DualComplex: The norm of this dual complex number
        """
        return DualComplex(self.a * self.a.conjugate(), self.a * self.b.conjugate() + self.b * self.a.conjugate())

    def inverse(self: 'DualComplex') -> 'DualComplex':
        """
        Calculates the inverse of this dual complex number.

            a^-1 - (a^-1 * b * a^-1) * ϵ
        
        Returns:
            DualComplex: The inverse of this dual complex number
        """
        a_inv = 1 / self.a
        return DualComplex(a_inv, -a_inv() * self.b * a_inv)

    @staticmethod
    def exp(Z: 'DualComplex') -> 'DualComplex':
        """
        Calculates e to the power of a dual complex number.

            e^(a + b * ϵ) = e^a * (1 + b * ϵ)
        
        Parameters:
            Z (DualComplex): The dual complex number
        
        Returns:
            DualComplex: e raised to the power of the dual complex number
        """
        return DualComplex(cm.exp(Z.a), cm.exp(Z.a) * Z.b)

    @staticmethod
    def log(Z: 'DualComplex') -> 'DualComplex':
        """
        Calculates the logarithm of a dual complex number.
        
            log(a + b * ϵ) = log(a) + (b / a) * ϵ

        Which is only defined when a != 0.
        
        Parameters:
            Z (DualComplex): The dual complex number

        Returns:
            DualComplex: The logarithm of the dual complex number
        """
        if Z.a == 0:
            raise ValueError("Log undefined: requires a != 0")
        else:
            return DualComplex(cm.log(Z.a), Z.b / Z.a)
        
    @staticmethod
    def sqrt(Z: 'DualComplex') -> 'DualComplex':
        """
        Calculates the square root of a dual complex number.

            sqrt(a + b * ϵ) = sqrt(a) + (b / (2 * sqrt(a))) * ϵ
        
        Which is only defined when a != 0.

        Parameters:
            Z (DualComplex): The dual complex number
        
        Returns:
            DualComplex: The square root of the dual complex number
        """
        if Z.a <= 0:
            raise ValueError("Sqrt undefined: requires a != 0")
        else:
            return DualComplex(cm.sqrt(Z.a), Z.b / (2 * cm.sqrt(Z.a)))
    
    @staticmethod
    def root(n: int, Z: 'DualComplex') -> 'DualComplex':
        """
        Calculates the n-th root of a dual complex number.

            root(n, Z) = exp((1 / n) * log(Z))
        
        Which is only defined when Z = 0 or a != 0.

        Expects n > 0.
        
        Parameters:
            Z (DualComplex): The dual complex number
        
        Returns:
            DualComplex: The n-th root of the dual complex number
        """
        if Z == 0:
            return DualComplex(0, 0)

        if n <= 0:
            raise ValueError("Root: expected n > 0")
        elif n == 1:
            return Z
        else:
            if Z.a == 0:
                raise ValueError("Root undefined: requires a != 0")
            else:
                return DualComplex.exp((1 / n) * DualComplex.log(Z))

class DualSplitComplex:
    """
    A class representing a dual split-complex number.
    
    A dual split-complex number is a number of the form
        a + b * ϵ
    with a, b in C' and ϵ^2 = 0.

    Attributes:
        a (SplitComplex): real component
        b (SplitComplex): dual component
    """

    def __init__(self: 'DualSplitComplex', a: SplitComplex, b: SplitComplex) -> 'DualSplitComplex':
        """
        Initializes a new dual split-complex number with the given components.

        Parameters:
            a (SplitComplex): real component
            b (SplitComplex): dual component
        
        Returns:
            DualSplitComplex: The new dual split-complex number
        """
        self.a = a
        self.b = b
    
    def __str__(self: 'DualSplitComplex') -> str:
        return f"({self.a}, {self.b})"

    def __repr__(self: 'DualSplitComplex') -> str:
        return f"DualSplitComplex({self.a}, {self.b})"
    
    def __eq__(self: 'DualSplitComplex', rhs) -> bool:
        if isinstance(rhs, (int, float)):
            return (self.a.a == rhs, self.a.b == 0 and self.b == 0)
        elif isinstance(rhs, SplitComplex):
            return (self.a == rhs and self.b == 0)
        elif isinstance(rhs, DualSplitComplex):
            return (self.a == rhs.a and self.b == rhs.b)
        else:
            return NotImplemented
    
    def __ne__(self: 'DualSplitComplex', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq

    def __abs__(self: 'DualSplitComplex') -> 'DualSplitComplex':
        return self.norm()
    
    def __neg__(self: 'DualSplitComplex') -> 'DualSplitComplex':
        return DualSplitComplex(-self.a, -self.b)
    
    def __add__(self: 'DualSplitComplex', rhs) -> 'DualSplitComplex':
        if isinstance(rhs, (int, float, SplitComplex)):
            return DualSplitComplex(self.a + rhs, self.b)
        elif isinstance(rhs, DualSplitComplex):
            return DualSplitComplex(self.a + rhs.a, self.b + rhs.b)
        else:
            return NotImplemented
    
    def __radd__(self: 'DualSplitComplex', lhs) -> 'DualSplitComplex':
        if isinstance(lhs, (int, float, SplitComplex)):
            return DualSplitComplex(lhs + self.a, self.b)
        elif isinstance(lhs, DualSplitComplex):
            return DualSplitComplex(lhs.a + self.a, lhs.b + self.b)
        else:
            return NotImplemented
    
    def __sub__(self: 'DualSplitComplex', rhs) -> 'DualSplitComplex':
        if isinstance(rhs, (int, float, SplitComplex)):
            return DualSplitComplex(self.a - rhs, self.b)
        elif isinstance(rhs, DualSplitComplex):
            return DualSplitComplex(self.a - rhs.a, self.b - rhs.b)
        else:
            return NotImplemented
    
    def __rsub__(self: 'DualSplitComplex', lhs) -> 'DualSplitComplex':
        if isinstance(lhs, (int, float, SplitComplex)):
            return DualSplitComplex(lhs - self.a, -self.b)
        elif isinstance(lhs, DualSplitComplex):
            return DualSplitComplex(lhs.a - self.a, lhs.b - self.b)
        else:
            return NotImplemented
    
    def __mul__(self: 'DualSplitComplex', rhs) -> 'DualSplitComplex':
        if isinstance(rhs, (int, float, SplitComplex)):
            return DualSplitComplex(self.a * rhs, self.b * rhs)
        elif isinstance(rhs, DualSplitComplex):
            return DualSplitComplex(self.a * rhs.a, self.a * rhs.b + self.b * rhs.a)
        else:
            return NotImplemented
    
    def __rmul__(self: 'DualSplitComplex', lhs) -> 'DualSplitComplex':
        if isinstance(lhs, (int, float, SplitComplex)):
            return DualSplitComplex(lhs * self.a, lhs * self.b)
        elif isinstance(lhs, DualSplitComplex):
            return DualSplitComplex(lhs.a * self.a, lhs.a * self.b + lhs.b * self.a)
        else:
            return NotImplemented
    
    def __truediv__(self: 'DualSplitComplex', rhs) -> 'DualSplitComplex':
        if isinstance(rhs, (int, float, SplitComplex)):
            return DualSplitComplex(self.a / rhs, self.b / rhs)
        elif isinstance(rhs, DualSplitComplex):
            return DualSplitComplex(self.a / rhs.a, (self.b * rhs.a - self.a * rhs.b) / rhs.a**2)
        else:
            return NotImplemented
    
    def __rtruediv__(self: 'DualSplitComplex', lhs) -> 'DualSplitComplex':
        if isinstance(lhs, (int, float, SplitComplex)):
            return DualSplitComplex(lhs / self.a, (-lhs * self.b) / self.a**2)
        elif isinstance(lhs, DualSplitComplex):
            return DualSplitComplex(lhs.a / self.a, (lhs.b * self.a - lhs.a * self.b) / self.a**2)
        else:
            return NotImplemented

    def __pow__(self: 'DualSplitComplex', exp) -> 'DualSplitComplex':
        if isinstance(exp, int):
            z = DualSplitComplex(SplitComplex(1, 0), SplitComplex(0, 0))

            if exp == 0:
                return z
            
            mul = self if exp > 0 else self.inverse()

            for _ in range(abs(exp)):
                z *= mul
            return z
        elif isinstance(exp, DualSplitComplex):
            if self.a.a <= abs(self.a.b):
                raise ValueError("Pow undefined: base dual split-complex number requires a > |b|")
            else:
                return DualSplitComplex.exp(exp * DualSplitComplex.log(self))
        else:
            return NotImplemented
    
    def __rpow__(self: 'DualSplitComplex', base) -> 'DualSplitComplex':
        if isinstance(base, (int, float)):
            if base <= 0:
                raise ValueError("Pow undefined: base requires x > 0")
            else:
                Z = DualSplitComplex(SplitComplex(base, 0), SplitComplex(0, 0))
                return Z**self
        elif isinstance(base, SplitComplex):
            if base.a <= abs(base.b):
                raise ValueError("Pow undefined: base split-complex number requires a > |b|")
            else:
                Z = DualSplitComplex(base, SplitComplex(0, 0))
                return Z**self
        elif isinstance(base, DualSplitComplex):
            return base**self
        else:
            return NotImplemented

    def conj_splitcomplex(self: 'DualSplitComplex') -> 'DualSplitComplex':
        """
        Returns the split-complex conjugate of this dual split-complex number.

        The split-complex conjugate of a dual split-complex number
            a + b * ϵ
        is
            a* + b* * ϵ

        Returns:
            DualSplitComplex: The split-complex conjugate of this dual split-complex number
        """
        return DualSplitComplex(self.a.conj(), self.b.conj())
    
    def conj_dual(self: 'DualSplitComplex') -> 'DualSplitComplex':
        """
        Returns the dual conjugate of this dual split-complex number.

        The dual conjugate of a dual split-complex number
            a + b * ϵ
        is
            a - b * ϵ

        Returns:
            DualSplitComplex: The dual conjugate of this dual split-complex number
        """
        return DualSplitComplex(self.a, -self.b)
    
    def conj_both(self: 'DualSplitComplex') -> 'DualSplitComplex':
        """
        Returns the full (split-complex and dual) conjugate of this dual split-complex number.

        The full conjugate of a dual split-complex number
            a + b * ϵ
        is
            a* - b* * ϵ

        Returns:
            DualSplitComplex: The full conjugate of this dual split-complex number
        """
        return DualSplitComplex(self.a.conj(), -self.b.conj())

    def norm(self: 'DualSplitComplex') -> 'DualSplitComplex':
        """
        Calculates the norm of this dual split-complex number.
        
            (a * a*) + (a * b* + b * a*) * ϵ

        Returns:
            DualSplitComplex: The norm of this dual split-complex number
        """
        return DualSplitComplex(self.a * self.a.conj(), self.a * self.b.conj() + self.b * self.a.conj())
    
    def inverse(self: 'DualSplitComplex') -> 'DualSplitComplex':
        """
        Calculates the inverse of this dual split-complex number.

            a^-1 - (a^-1 * b * a^-1) * ϵ
        
        Returns:
            DualSplitComplex: The inverse of this dual split-complex number
        """
        return DualSplitComplex(self.a.inverse(), -self.a.inverse() * self.b * self.a.inverse())
    
    @staticmethod
    def exp(Z: 'DualSplitComplex') -> 'DualSplitComplex':
        """
        Calculates e to the power of a dual split-complex number.

            e^(a + b * ϵ) = e^a * (1 + b * ϵ)
        
        Parameters:
            Z (DualSplitComplex): The dual split-complex number
        
        Returns:
            DualSplitComplex: e raised to the power of the dual split-complex number
        """
        return DualSplitComplex(SplitComplex.exp(Z.a), SplitComplex.exp(Z.a) * Z.b)

    @staticmethod
    def log(Z: 'DualSplitComplex') -> 'DualSplitComplex':
        """
        Calculates the logarithm of a dual split-complex number.
        
            log(a + b * ϵ) = log(a) + (b / a) * ϵ

        Which is only defined when Re(a) > |Im(a)|.
        
        Parameters:
            Z (DualSplitComplex): The dual split-complex number

        Returns:
            DualSplitComplex: The logarithm of the dual split-complex number
        """
        if Z.a.a <= abs(Z.a.b):
            raise ValueError("Log undefined: requires Re(a) > |Im(a)|")
        else:
            return DualSplitComplex(SplitComplex.log(Z.a), Z.b / Z.a)
    
    @staticmethod
    def sqrt(Z: 'DualSplitComplex') -> 'DualSplitComplex':
        """
        Calculates the square root of a dual split-complex number.

            sqrt(a + b * ϵ) = sqrt(a) + (b / (2 * sqrt(a))) * ϵ
        
        Which is only defined when Re(a) > |Im(a)|.

        Parameters:
            Z (DualSplitComplex): The dual split-complex number
        
        Returns:
            DualSplitComplex: The square root of the dual split-complex number
        """
        if Z.a.a <= abs(Z.a.b):
            raise ValueError("Sqrt undefined: requires Re(a) > |Im(a)|")
        else:
            return DualSplitComplex(SplitComplex.sqrt(Z.a), Z.b / (2 * SplitComplex.sqrt(Z.a)))
    
    @staticmethod
    def root(n: int, Z: 'DualSplitComplex') -> 'DualSplitComplex':
        """
        Calculates the n-th root of a dual split-complex number.

            root(n, Z) = exp((1 / n) * log(Z))
        
        Which is only defined when Z = 0 or Re(a) > |Im(a)|.

        Expects n > 0.
        
        Parameters:
            Z (DualSplitComplex): The dual split-complex number
        
        Returns:
            DualSplitComplex: The n-th root of the dual split-complex number
        """
        if Z == 0:
            return DualSplitComplex(0, 0)

        if n <= 0:
            raise ValueError("Root: expected n > 0")
        elif n == 1:
            return Z
        else:
            if Z.a.a <= abs(Z.a.b):
                raise ValueError("Root undefined: requires Re(a) > |Im(a)|")
            else:
                return DualSplitComplex.exp((1 / n) * DualSplitComplex.log(Z))

class DualQuat:
    """
    A class representing a dual quaternion.
    
    A dual quaternion is a number of the form
        a + b * ϵ
    with a, b in H and ϵ^2 = 0.

    Attributes:
        a (Quat): real component
        b (Quat): dual component
    """

    def __init__(self: 'DualQuat', a: Quat, b: Quat) -> 'DualQuat':
        """
        Initializes a new dual quaternion with the given components.

        Parameters:
            a (Quat): real component
            b (Quat): dual component
        
        Returns:
            DualQuat: The new dual quaternion
        """
        self.a = a
        self.b = b
    
    def __str__(self: 'DualQuat') -> str:
        return f"({self.a}, {self.b})"

    def __repr__(self: 'DualQuat') -> str:
        return f"DualQuat({self.a}, {self.b})"
    
    def __eq__(self: 'DualQuat', rhs) -> bool:
        if isinstance(rhs, (int, float)):
            return (self.a.w == rhs and self.a.x == 0 and self.a.y == 0 and self.a.z == 0 and self.b == 0)
        elif isinstance(rhs, Quat):
            return (self.a == rhs and self.b == 0)
        elif isinstance(rhs, DualQuat):
            return (self.a == rhs.a and self.b == rhs.b)
        else:
            return NotImplemented

    def __ne__(self: 'DualQuat', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq
    
    def __abs__(self: 'DualQuat') -> 'DualQuat':
        return self.norm()

    def __neg__(self: 'DualQuat') -> 'DualQuat':
        return DualQuat(-self.a, -self.b)

    def __add__(self: 'DualQuat', rhs) -> 'DualQuat':
        if isinstance(rhs, (int, float, Quat)):
            return DualQuat(self.a + rhs, self.b)
        elif isinstance(rhs, DualQuat):
            return DualQuat(self.a + rhs.a, self.b + rhs.b)
        else:
            return NotImplemented
    
    def __radd__(self: 'DualQuat', lhs) -> 'DualQuat':
        if isinstance(lhs, (int, float, Quat)):
            return DualQuat(lhs + self.a, self.b)
        elif isinstance(lhs, DualQuat):
            return DualQuat(lhs.a + self.a, lhs.b + self.b)
        else:
            return NotImplemented
    
    def __sub__(self: 'DualQuat', rhs) -> 'DualQuat':
        if isinstance(rhs, (int, float, Quat)):
            return DualQuat(self.a - rhs, self.b)
        elif isinstance(rhs, DualQuat):
            return DualQuat(self.a - rhs.a, self.b - rhs.b)
        else:
            return NotImplemented
    
    def __rsub__(self: 'DualQuat', lhs) -> 'DualQuat':
        if isinstance(lhs, (int, float, Quat)):
            return DualQuat(lhs - self.a, -self.b)
        elif isinstance(lhs, DualQuat):
            return DualQuat(lhs.a - self.a, lhs.b - self.b)
        else:
            return NotImplemented
    
    def __mul__(self: 'DualQuat', rhs) -> 'DualQuat':
        if isinstance(rhs, (int, float, Quat)):
            return DualQuat(self.a * rhs, self.b * rhs)
        elif isinstance(rhs, DualQuat):
            return DualQuat(self.a * rhs.a, self.a * rhs.b + self.b * rhs.a)
        else:
            return NotImplemented
    
    def __rmul__(self: 'DualQuat', lhs) -> 'DualQuat':
        if isinstance(lhs, (int, float, Quat)):
            return DualQuat(lhs * self.a, lhs * self.b)
        elif isinstance(lhs, DualQuat):
            return DualQuat(lhs.a * self.a, lhs.a * self.b + lhs.b * self.a)
        else:
            return NotImplemented
        
    # NOTE: __truediv__ and __rtruediv__ are emitted to force multiplication with the inverse which should aid with clarity

    def __pow__(self: 'DualQuat', exp) -> 'DualQuat':
        if isinstance(exp, int):
            q = DualQuat(Quat(1, 0, 0, 0), Quat(0, 0, 0, 0))

            if exp == 0:
                return q
            
            mul = self if exp > 0 else self.inverse()

            for _ in range(abs(exp)):
                q *= mul
            return q
        elif isinstance(exp, DualQuat):
            if self.a == 0:
                raise ValueError("Pow undefined: base dual quaternion requires a != 0")
            else:
                return DualQuat.exp(exp * DualQuat.log(self))
        else:
            return NotImplemented

    def __rpow__(self: 'DualQuat', base) -> 'DualQuat':
        if isinstance(base, (int, float)):
            if base == 0:
                raise ValueError("Pow undefined: base requires x != 0")
            else:
                Q = DualQuat(Quat(base, 0, 0, 0), Quat(0, 0, 0, 0))
                return Q**self
        elif isinstance(base, Quat):
            if base == 0:
                raise ValueError("Pow undefined: base quaternion requires q != 0")
            else:
                Q = DualQuat(base, Quat(0, 0, 0, 0))
                return Q**self
        elif isinstance(base, DualQuat):
            return base**self
        else:
            return NotImplemented

    def conj_quat(self: 'DualQuat') -> 'DualQuat':
        """
        Returns the quaternion conjugate of this dual quaternion.

        The quaternion conjugate of a dual quaternion
            a + b * ϵ
        is
            a* + b* * ϵ

        Returns:
            DualQuat: The quaternion conjugate of this dual quaternion
        """
        return DualQuat(self.a.conj(), self.b.conj())

    def conj_dual(self: 'DualQuat') -> 'DualQuat':
        """
        Returns the dual conjugate of this dual quaternion.

        The dual conjugate of a dual quaternion
            a + b * ϵ
        is
            a - b * ϵ

        Returns:
            DualQuat: The dual conjugate of this dual quaternion
        """
        return DualQuat(self.a, -self.b)

    def conj_both(self: 'DualQuat') -> 'DualQuat':
        """
        Returns the full (quaternion and dual) conjugate of this dual quaternion.

        The full conjugate of a dual quaternion
            a + b * ϵ
        is
            a* - b* * ϵ

        Returns:
            DualQuat: The full conjugate of this dual quaternion
        """
        return DualQuat(self.a.conj(), -self.b.conj())
    
    def norm(self: 'DualQuat') -> 'DualQuat':
        """
        Calculates the norm of this dual quaternion.
        
            (a * a*) + (a * b* + b * a*) * ϵ

        Returns:
            DualQuat: The norm of this dual quaternion
        """
        return DualQuat(self.a * self.a.conj(), self.a * self.b.conj() + self.b * self.a.conj())

    def inverse(self: 'DualQuat') -> 'DualQuat':
        """
        Calculates the inverse of this dual quaternion.

            a^-1 - (a^-1 * b * a^-1) * ϵ
        
        Returns:
            DualQuat: The inverse of this dual quaternion
        """
        return DualQuat(self.a.inverse(), -self.a.inverse() * self.b * self.a.inverse())

    @staticmethod
    def from_vec3(vec: Vec3) -> 'DualQuat':
        """
        Initializes a dual quaternion from a Vec3 for use in rotation and translation.

        A Vec3
            (v_x, v_y, v_z)
        results in
            1 + (v_x * i + v_y * j + v_z * k) ϵ

        Parameters:
            vec (Vec3): The 3-dimensional vector
        
        Returns:
            DualQuat: A dual quaternion that represents the Vec3
        """
        return DualQuat(Quat(1, 0, 0, 0), Quat.from_vec3(vec))

    @staticmethod
    def exp(Q: 'DualQuat') -> 'DualQuat':
        """
        Calculates e to the power of a dual quaternion.

            e^(a + b * ϵ) = e^a * (1 + b * ϵ)
        
        Parameters:
            Q (DualQuat): The dual quaternion
        
        Returns:
            DualQuat: e raised to the power of the dual quaternion
        """
        return DualQuat(Quat.exp(Q.a), Quat.exp(Q.a) * Q.b)
    
    @staticmethod
    def log(Q: 'DualQuat') -> 'DualQuat':
        """
        Calculates the logarithm of a dual quaternion.
        
            log(a + b * ϵ) = log(a) + (b / a) * ϵ

        Which is only defined when a != 0.
        
        Parameters:
            Q (DualQuat): The dual quaternion

        Returns:
            DualQuat: The logarithm of the dual quaternion
        """
        if Q.a == 0:
            raise ValueError("Log undefined: requires a != 0")
        else:
            return DualQuat(Quat.log(Q.a), Q.b * Q.a.inverse())
    
    @staticmethod
    def sqrt(Q: 'DualQuat') -> 'DualQuat':
        """
        Calculates the square root of a dual quaternion.

            sqrt(a + b * ϵ) = sqrt(a) + (b / (2 * sqrt(a))) * ϵ
        
        Which is only defined when a != 0.

        Parameters:
            Q (DualQuat): The dual quaternion
        
        Returns:
            DualQuat: The square root of the dual quaternion
        """
        if Q.a == 0:
            raise ValueError("Sqrt undefined: requires a != 0")
        else:
            return DualQuat(Quat.sqrt(Q.a), Q.b * (.5 * Quat.sqrt(Q.a).inverse()))

    @staticmethod
    def root(n: int, Q: 'DualQuat') -> 'DualQuat':
        """
        Calculates the n-th root of a dual quaternion.

            root(n, Q) = exp((1 / n) * log(Q))
        
        Which is only defined when Q = 0 or a != 0.

        Expects n > 0.
        
        Parameters:
            Q (DualQuat): The dual quaternion
        
        Returns:
            DualQuat: The n-th root of the dual quaternion
        """
        if Q == 0:
            return DualQuat(0, 0)

        if n <= 0:
            raise ValueError("Root: expected n > 0")
        elif n == 1:
            return Q
        else:
            if Q.a == 0:
                raise ValueError("Root undefined: requires a != 0")
            else:
                return DualQuat.exp((1 / n) * DualQuat.log(Q))

    @staticmethod
    def rot_trans(rot: Quat, trans: Quat) -> 'DualQuat':
        """
        Helper function that creates a dual quaternion to be used for rotation and translation.

        Parameters:
            rot (Quat): The rotation quaternion with length 1
            trans (Quat): The imaginary translation quaternion
        
        Returns:
            DualQuat: The dual quaternion that applies the rotation and translation
        """
        return DualQuat(rot, .5 * trans * rot)

    @staticmethod
    def rot_trans2(vec: Vec3, phi: float, trans: Quat):
        """
        Helper function that creates a dual quaternion to be used for rotation and translation.

        Parameters:
            vec (Vec3): The rotation axis with length 1 (normalized)
            phi (float): The rotation angle
            trans (Quat): The imaginary translation quaternion
        
        Returns:
            DualQuat: The dual quaternion that applies the rotation and translation
        """
        rot = Quat.exp_rot(vec, phi)
        return DualQuat.rot_trans(rot, trans)

class DualSplitQuat:
    """
    A class representing a dual split-quaternion.
    
    A dual split-quaternion is a number of the form
        a + b * ϵ
    with a, b in H' and ϵ^2 = 0.

    Attributes:
        a (SplitQuat): real component
        b (SplitQuat): dual component
    """

    def __init__(self: 'DualSplitQuat', a: SplitQuat, b: SplitQuat) -> 'DualSplitQuat':
        """
        Initializes a new dual split-quaternion with the given components.

        Parameters:
            a (SplitQuat): real component
            b (SplitQuat): dual component
        
        Returns:
            DualSplitQuat: The new dual split-quaternion
        """
        self.a = a
        self.b = b
    
    def __str__(self: 'DualSplitQuat') -> str:
        return f"({self.a}, {self.b})"
    
    def __repr__(self: 'DualSplitQuat') -> str:
        return f"DualSplitQuat({self.a}, {self.b})"
    
    def __eq__(self: 'DualSplitQuat', rhs) -> bool:
        if isinstance(rhs, (int, float)):
            return (self.a.w == rhs and self.a.x == 0 and self.a.y == 0 and self.a.z == 0 and self.b == 0)
        elif isinstance(rhs, SplitQuat):
            return (self.a == rhs and self.b == 0)
        elif isinstance(rhs, DualSplitQuat):
            return (self.a == rhs.a and self.b == rhs.b)
        else:
            return NotImplemented

    def __ne__(self: 'DualSplitQuat', rhs) -> bool:
        eq = self.__eq__(rhs)
        return NotImplemented if eq is NotImplemented else not eq

    def __abs__(self: 'DualSplitQuat') -> 'DualSplitQuat':
        return self.norm()
    
    def __neg__(self: 'DualSplitQuat') -> 'DualSplitQuat':
        return DualSplitQuat(-self.a, -self.b)
    
    def __add__(self: 'DualSplitQuat', rhs) -> 'DualSplitQuat':
        if isinstance(rhs, (int, float, SplitQuat)):
            return DualSplitQuat(self.a + rhs, self.b)
        elif isinstance(rhs, DualSplitQuat):
            return DualSplitQuat(self.a + rhs.a, self.b + rhs.b)
        else:
            return NotImplemented
    
    def __radd__(self: 'DualSplitQuat', lhs) -> 'DualSplitQuat':
        if isinstance(lhs, (int, float, SplitQuat)):
            return DualSplitQuat(lhs + self.a, self.b)
        elif isinstance(lhs, DualSplitQuat):
            return DualSplitQuat(lhs.a + self.a, lhs.b + self.b)
        else:
            return NotImplemented
    
    def __sub__(self: 'DualSplitQuat', rhs) -> 'DualSplitQuat':
        if isinstance(rhs, (int, float, SplitQuat)):
            return DualSplitQuat(self.a - rhs, self.b)
        elif isinstance(rhs, DualSplitQuat):
            return DualSplitQuat(self.a - rhs.a, self.b - rhs.b)
        else:
            return NotImplemented
    
    def __rsub__(self: 'DualSplitQuat', lhs) -> 'DualSplitQuat':
        if isinstance(lhs, (int, float, SplitQuat)):
            return DualSplitQuat(lhs - self.a, -self.b)
        elif isinstance(lhs, DualSplitQuat):
            return DualSplitQuat(lhs.a - self.a, lhs.b - self.b)
        else:
            return NotImplemented

    def __mul__(self: 'DualSplitQuat', rhs) -> 'DualSplitQuat':
        if isinstance(rhs, (int, float, SplitQuat)):
            return DualSplitQuat(self.a * rhs, self.b * rhs)
        elif isinstance(rhs, DualSplitQuat):
            return DualSplitQuat(self.a * rhs.a, self.a * rhs.b + self.b * rhs.a)
        else:
            return NotImplemented
    
    def __rmul__(self: 'DualSplitQuat', lhs) -> 'DualSplitQuat':
        if isinstance(lhs, (int, float, SplitQuat)):
            return DualSplitQuat(lhs * self.a, lhs * self.b)
        elif isinstance(lhs, DualSplitQuat):
            return DualSplitQuat(lhs.a * self.a, lhs.a * self.b + lhs.b * self.a)
        else:
            return NotImplemented
    
    # NOTE: __truediv__ and __rtruediv__ are emitted to force multiplication with the inverse which should aid with clarity

    def __pow__(self: 'DualSplitQuat', exp) -> 'DualSplitQuat':
        if isinstance(exp, int):
            q = DualSplitQuat(SplitQuat(1, 0, 0, 0), SplitQuat(0, 0, 0, 0))

            if exp == 0:
                return q
            
            mul = self if exp > 0 else self.inverse()

            for _ in range(abs(exp)):
                q *= mul
            return q
        elif isinstance(exp, DualSplitQuat):
            if self.a == 0 or self.a.norm() <= 0:
                raise ValueError("Pow undefined: base dual split-quaternion requires a != 0 and |a| > 0")
            else:
                return DualSplitQuat.exp(exp * DualSplitQuat.log(self))
        else:
            return NotImplemented
    
    def __rpow__(self: 'DualSplitQuat', base) -> 'DualSplitQuat':
        if isinstance(base, (int, float)):
            if base <= 0:
                raise ValueError("Pow undefined: base requires x > 0")
            else:
                Q = DualSplitQuat(SplitQuat(base, 0, 0, 0), SplitQuat(0, 0, 0, 0))
                return Q**self
        elif isinstance(base, SplitQuat):
            if base == 0 or base.norm() <= 0:
                raise ValueError("Pow undefined: base split-quaternion requires q != 0 and |q| > 0")
            else:
                Q = DualSplitQuat(base, SplitQuat(0, 0, 0, 0))
                return Q**self
        elif isinstance(base, DualSplitQuat):
            return base**self
        else:
            return NotImplemented

    def conj_splitquat(self: 'DualSplitQuat') -> 'DualSplitQuat':
        """
        Returns the split-quaternion conjugate of this dual split-quaternion.

        The split-quaternion conjugate of a dual split-quaternion
            a + b * ϵ
        is
            a* + b* * ϵ

        Returns:
            DualSplitQuat: The split-quaternion conjugate of this dual split-quaternion
        """
        return DualSplitQuat(self.a.conj(), self.b.conj())
    
    def conj_dual(self: 'DualSplitQuat') -> 'DualSplitQuat':
        """
        Returns the dual conjugate of this dual split-quaternion.

        The dual conjugate of a dual split-quaternion
            a + b * ϵ
        is
            a - b * ϵ

        Returns:
            DualSplitQuat: The dual conjugate of this dual split-quaternion
        """
        return DualSplitQuat(self.a, -self.b)
    
    def conj_both(self: 'DualSplitQuat') -> 'DualSplitQuat':
        """
        Returns the full (split-quaternion and dual) conjugate of this dual split-quaternion.

        The full conjugate of a dual split-quaternion
            a + b * ϵ
        is
            a* - b* * ϵ

        Returns:
            DualSplitQuat: The full conjugate of this dual split-quaternion
        """
        return DualSplitQuat(self.a.conj(), -self.b.conj())

    def norm(self: 'DualSplitQuat') -> 'DualSplitQuat':
        """
        Calculates the norm of this dual split-quaternion.
        
            (a * a*) + (a * b* + b * a*) * ϵ

        Returns:
            DualSplitQuat: The norm of this dual split-quaternion
        """
        return DualSplitQuat(self.a * self.a.conj(), self.a * self.b.conj() + self.b * self.a.conj())

    def inverse(self: 'DualSplitQuat') -> 'DualSplitQuat':
        """
        Calculates the inverse of this dual split-quaternion.

            a^-1 - (a^-1 * b * a^-1) * ϵ
        
        Returns:
            DualSplitQuat: The inverse of this dual split-quaternion
        """
        return DualSplitQuat(self.a.inverse(), -self.a.inverse() * self.b * self.a.inverse())
    
    @staticmethod
    def exp(Q: 'DualSplitQuat') -> 'DualSplitQuat':
        """
        Calculates e to the power of a dual split-quaternion.

            e^(a + b * ϵ) = e^a * (1 + b * ϵ)
        
        Parameters:
            Q (DualSplitQuat): The dual split-quaternion
        
        Returns:
            DualSplitQuat: e raised to the power of the dual split-quaternion
        """
        return DualSplitQuat(SplitQuat.exp(Q.a), SplitQuat.exp(Q.a) * Q.b)

    @staticmethod
    def log(Q: 'DualSplitQuat') -> 'DualSplitQuat':
        """
        Calculates the logarithm of a dual split-quaternion.
        
            log(a + b * ϵ) = log(a) + (b / a) * ϵ

        Which is only defined when a != 0 and |a| > 0.
        
        Parameters:
            Q (DualSplitQuat): The dual split-quaternion

        Returns:
            DualSplitQuat: The logarithm of the dual split-quaternion
        """
        if Q.a == 0 or Q.a.norm() <= 0:
            raise ValueError("Log undefined: requires a != 0 and |a| > 0")
        else:
            return DualSplitQuat(SplitQuat.log(Q.a), Q.b * Q.a.inverse())
    
    @staticmethod
    def sqrt(Q: 'DualSplitQuat') -> 'DualSplitQuat':
        """
        Calculates the square root of a dual split-quaternion.

            sqrt(a + b * ϵ) = sqrt(a) + (b / (2 * sqrt(a))) * ϵ
        
        Which is only defined when a != 0 and |a| > 0.

        Parameters:
            Q (DualSplitQuat): The dual split-quaternion
        
        Returns:
            DualSplitQuat: The square root of the dual split-quaternion
        """
        if Q.a == 0 or Q.a.norm() <= 0:
            raise ValueError("Sqrt undefined: requires a != 0 and |a| > 0")
        else:
            return DualSplitQuat(SplitQuat.sqrt(Q.a), Q.b * (.5 * SplitQuat.sqrt(Q.a).inverse()))

    @staticmethod
    def root(n: int, Q: 'DualSplitQuat') -> 'DualSplitQuat':
        """
        Calculates the n-th root of a dual split-quaternion.

            root(n, Q) = exp((1 / n) * log(Q))
        
        Which is only defined when Q = 0 or a != 0 and |a| > 0.

        Expects n > 0.
        
        Parameters:
            Q (DualSplitQuat): The dual split-quaternion
        
        Returns:
            DualSplitQuat: The n-th root of the dual split-quaternion
        """
        if Q == 0:
            return DualSplitQuat(0, 0)

        if n <= 0:
            raise ValueError("Root: expected n > 0")
        elif n == 1:
            return Q
        else:
            if Q.a == 0 or Q.a.norm() <= 0:
                raise ValueError("Root undefined: requires a != 0 and |a| > 0")
            else:
                return DualSplitQuat.exp((1 / n) * DualSplitQuat.log(Q))

    @staticmethod
    def boost_transI(rapidity: float, trans: SplitQuat) -> 'DualSplitQuat':
        """
        Helper function that creates a dual split-quaternion to be used for a boost in the i-direction and translation.

        Parameters:
            rapidity (float): The rapidity
            trans (SplitQuat): The imaginary translation split-quaternion
        
        Returns:
            DualSplitQuat: The dual split-quaternion that applies the boost in the i-direction and translation
        """
        boost = SplitQuat.expI(rapidity)
        return DualSplitQuat(boost, .5 * trans * boost)

    @staticmethod
    def boost_transJ(rapidity: float, trans: SplitQuat) -> 'DualSplitQuat':
        """
        Helper function that creates a dual split-quaternion to be used for a boost in the j-direction and translation.

        Parameters:
            rapidity (float): The rapidity
            trans (SplitQuat): The imaginary translation split-quaternion
        
        Returns:
            DualSplitQuat: The dual split-quaternion that applies the boost in the j-direction and translation
        """
        boost = SplitQuat.expJ(rapidity)
        return DualSplitQuat(boost, .5 * trans * boost)
    
    @staticmethod
    def rot_transK(phi: float, trans: SplitQuat) -> 'DualSplitQuat':
        """
        Helper function that creates a dual split-quaternion to be used for a rotation in the k-direction and translation.

        Parameters:
            phi (float): The rotation angle
            trans (SplitQuat): The imaginary translation split-quaternion
        
        Returns:
            DualSplitQuat: The dual split-quaternion that applies the rotation in the k-direction and translation
        """
        rot = SplitQuat.expK(phi)
        return DualSplitQuat(rot, .5 * trans * rot)

    @staticmethod
    def lorentz_trans(lor: SplitQuat, trans: SplitQuat) -> 'DualSplitQuat':
        """
        Helper function that creates a dual split-quaternion representing a general Lorentz transformation and translation.

        Parameters:
            lor (SplitQuat): The Lorentz transformation split-quaternion with length 1
            trans (SplitQuat): The imaginary translation split-quaternion
        
        Returns:
            DualSplitQuat: The dual split-quaternion that applies the Lorentz transformation and translation
        """
        return DualSplitQuat(lor, .5 * trans * lor)