import math

class MyMatrix:

    def __init__(self, A):
        """
        To be consistent with code already written, it is
        being assumed that the matrix A (which was used
        because it is a common name for a 'general'  matrix
        in linear algebra) is a 2-D array, or nested list
        """
        self.A = A
  
    def __str__(self):
        pass
  
    def __mul__(self, A, B):
        """
        Takes in two MyMatrix objects A and B.
        Returns AB as a MyMatrix object. I found it made the most
        sense to code this by transposing B and then simply dotting
        the rows of A and B^T, but if there's a more efficient
        way that would be ideal.
        """
        if len(A[0]) != len(B):
            return None

        ab = [[0 for j in range(len(B[0]))] for i in range(len(A))]
        bt = self.transpose(B)

        for i in range(len(A)):
            for j in range(len(bt)):
                ab[i][j] = self._dot(A[i], bt[j])

        return MyMatrix(ab)

    def __len__(self):
        return len(self.A)

    def __getitem__(self, row):
        return self.A[row]

    def transpose(self, A):
        """
        Takes in A, a MyMatrix object.
        Returns the transpose of A as a MyMatrix object.
        """
        at = [[0 for j in range(len(A))] for i in range(len(A[0]))]

        for i in range(len(A)):
            for j in range(len(A[0])):
                at[j][i] = A[i][j]

        return MyMatrix(at)

    def _dot(self, a, b):
        """
        Returns the dot product of two vectors (represented
        as 1-D arrays or lists). Used primarily in matrix
        multiplication in the MyMatrix class, and so it's been
        made 'private' as this isn't a MyVector class.
        """
        ab = 0
        for i in range(len(a)):
            ab += (a[i] * b[i])

        return ab

  
    def __add__(self, x, y):
        """
        x, y are MyMatrix object, the return is also MyMatrix object.
        here, I assume __getitem__ is defined
    	here, I assume the init of MyMatrix using the nested list,
    	which is like this [[], []]
        """
        sums = []
        for i in range(len(x)):
        	row = []
        	for j in range(len(y)):
        		# row.append(x.index(i, j) + y.index(i, j))
        		row.append(x[i][j] + y[i][j])
        	sums.append(row)
    
        return MyMatrix(sums)

    def tri_lower(self):
        """
        Returns the lower triangular portion of the matrix. All other values are
        zeroed out.
        """
        au = self
        k = self.transpose(au)
        k2 = len(k)
        r = 1
        for p in range(0, len(self)):
            for q in range(r, k2):
                au[p][q] = 0
            r += 1
        return au

    def norm1(self):
        """
        A is MyMatrix object which is a list in the form [[],[],[]]
        The method returns the 1-norm of matrix A, ||A||_1, which is a number
        The formula to calculate 1 norm of a matrix is ||A||_1 = max_{1<=j<=n} sum_{i=1^m} |a_{ij}|
        which is maximum absolute column sum of the matrix
        Example: input: A = MyMatrix([[1,2,3],[4,5,6]])
                 return: A.norm1() = 9
        """
        matrix_norm1 = 0
        
        for j in range(0,len(self.A[0])):#the number of columns
            column_sum = 0
            for i in range(0,len(self.A)):#the number of rows
               column_sum += abs(self.A[i][j])
            if matrix_norm1 < column_sum:
                matrix_norm1 = column_sum
            else:
                continue
        
        return matrix_norm1

    def norm_infinity(self):
        """
        A is MyMatrix object which is a list in the form [[],[],[]]
        The method returns the infinity-norm of matrix A, ||A||_∞, which is a number
        The formula to calculate infinity norm of a matrix is ||A||_∞ = max_{1<=i<=n} sum_{j=1^n} |a_{ij}|
        which is maximum absolute row sum of the matrix
        Example: input: A = MyMatrix([[1,2,3],[4,5,6]])
                 return: A.norm_infinity() = 15
        """
        
        matrix_norm_infinity = 0
        
        for i in range(0, len(self.A)):
            row_sum = 0
            for j in range(0, len(self.A[0])):
                row_sum += abs(self.A[i][j])
            if matrix_norm_infinity < row_sum:
                matrix_norm_infinity = row_sum
            else:
                continue
        
        return matrix_norm_infinity

    def normF(self):
        """
        A is MyMatrix object which is a list in the form [[],[],[]]
        The method returns the Frobenius norm of matrix A, ||A||_F, which is a number
        The formula to calculate the Frobenius norm of a matrix is ||A||_F = sqrt(sum_{1<=i<=n} (sum_{j=1^n} |a_{ij}|^2))
        Example: input: A = MyMatrix([[1,0,1],[1,2,0]])
                 return: A.normF() = sqrt(7)
        """
        
        matrix_normF = 0
        
        for i in range(0, len(self.A)):
            for j in range(0, len(self.A[0])):
                matrix_normF += abs(self.A[i][j])**2
        
        return math.sqrt(matrix_normF)
    
    def solveUpperTri(self, b):
        """
        This method returns the solution vector for an upper triangular system. 
        If the length of the matrix and vector are not correct, it will return none.
        """
        
        if ((len(b) != (len(self.A[0])))):
            return None;
    
        myNewVector = [0 for i in range(len(self.A))]
    
        for i in range(len(self.A)-1,-1, -1):
            mySum = b[i]
    
            for j in range(len(self.A)-1,i, -1):
                myDiff = (self.A[i][j]*myNewVector[j])
                mySum = mySum - myDiff
    
            myNewVector[i] = mySum/self.A[i][i]
    
        return myNewVector


if __name__  == "__main__":
    # some test examples
    A1 = MyMatrix([[1,2,3],[3,4,5]])
    B1 = MyMatrix([[1,2,3,4],[5,6,7,8],[9,10,11,12]])
    myA = MyMatrix([[1,2,3],[0,-1,-1],[0,0,-1]])
    myB = MyMatrix([1,1,-1])
    myX = myA.solveUpperTri(myB)
    print(myX)
    
    
    
 # Start: KC Gubler Edit
    
from sympy import Symbol, Derivative
import math
import sympy as sym

x = Symbol('x')

y = Symbol('y')

cos = Symbol('cos()')

sin = Symbol('sin()')

math.e = Symbol('e')

def jacobian(u,v):

   # u = 2*x*y**2 + 4*y**3 + x
   # v = 2*y**5 + x*y

    dxdu = Derivative(u, x)
    dydu = Derivative(u, y)
    dxdv = Derivative(v, x)
    dydv = Derivative(v, y)

    jacobian = dxdu*dydv - dxdv*dydu
    jacobian = jacobian.doit()

    return jacobian


jacobian(sym.sin(x), (2*y**5 + x*y))
jacobian(math.e**x, y)

# End: KC Gubler Edit
