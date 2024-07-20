import numpy as np
import sympy as sp

# Function to convert a matrix to a string representation
def matrix_to_string(matrix):
    matrix_str = "["
    for i, row in enumerate(matrix):
        row_str = " ".join(f"{int(x)}" if x.is_integer() else f"{x}" for x in row)
        matrix_str += f"[{row_str}]"
        if i < matrix.shape[0] - 1:
            matrix_str += "\n "
    matrix_str += "]"
    return matrix_str

# Function to fix values near zero
def fix_near_zero(value, tolerance=1e-4):
    if np.isclose(value, 0 , atol=tolerance):
        return 0.0
    elif np.isclose(value, 1, atol=tolerance):
        return 1.0
    return value

# Function to perform type II and III row operations to obtain RREF
def typeII_and_III_RREF(matrix):
    
    x = sp.symbols('x')
    pivot_indices = []
    numberOfRows, numberOfColumns = matrix.shape

    # The first non-zero entry of each row must be 1
    for i in range(numberOfRows):
        first_nonZero = False  # Used to flag if we found the pivot
        for jth_column, entry in enumerate(matrix[i, :]):
            if not first_nonZero:
                # If this is a non-zero element, divide the entire row by it to make it a leading 1 and mark that we found the pivot
                entry = fix_near_zero(entry)
                if entry != 0:
                    matrix[i, :] /= entry
                    first_nonZero = True
                    pivot_row_index = i 
                    pivot_column_index = jth_column
                    pivot_indices.append((pivot_row_index, pivot_column_index))

                    # Type III: Make the other entries in this column zero
                    for ith_row, column_entry in enumerate(matrix[:, jth_column]):
                        if (ith_row != i) and (column_entry != 0):
                            equation = sp.Eq(matrix[ith_row][jth_column] + x * matrix[i, jth_column], 0)
                            solution = sp.solve(equation, x)
                            float_solution = float(solution[0])
                            matrix[ith_row, :] += float_solution * matrix[i, :]

                else:
                    pass
            else:
                break

    return pivot_indices

# Function to perform type II and III row operations to obtain inverse
def typeII_and_III_Inverse(matrix, identity_matrix):
    
    x = sp.symbols('x')
    pivot_indices = []
    numberOfRows, numberOfColumns = matrix.shape

    # The first non-zero entry of each row must be 1
    for i in range(numberOfRows):
        first_nonZero = False  # Used to flag if we found the pivot
        for jth_column, entry in enumerate(matrix[i, :]):
            if not first_nonZero:
                # If this is a non-zero element, divide the entire row by it to make it a leading 1 and mark that we found the pivot
                entry = fix_near_zero(entry)
                if entry != 0:
                    matrix[i, :] /= entry
                    identity_matrix[i, :] /= entry
                    first_nonZero = True
                    pivot_row_index = i 
                    pivot_column_index = jth_column
                    pivot_indices.append((pivot_row_index, pivot_column_index))

                    # Type III: Make the other entries in this column zero
                    for ith_row, column_entry in enumerate(matrix[:, jth_column]):
                        if (ith_row != i) and (column_entry != 0):
                            equation = sp.Eq(matrix[ith_row][jth_column] + x * matrix[i, jth_column], 0)
                            solution = sp.solve(equation, x)
                            float_solution = float(solution[0])
                            matrix[ith_row, :] += float_solution * matrix[i, :]
                            identity_matrix[ith_row, :] += float_solution * identity_matrix[i, :]

                else:
                    pass
            else:
                break

    return pivot_indices

# Function to perform type I row operations to put pivots in order
def typeI_pivots_RREF(matrix, pivot_index):
    pivot_indices = pivot_index
    new_pivot_indices = pivot_indices
    while True:
        swapped = False
        pivot_indices = new_pivot_indices

        for stopped_pivot in pivot_indices:
            for moving_pivot in pivot_indices:
                if ((stopped_pivot[0] > moving_pivot[0]) and (stopped_pivot[1] < moving_pivot[1])) or \
                ((stopped_pivot[0] < moving_pivot[0]) and (stopped_pivot[1] > moving_pivot[1])) and ((stopped_pivot[0] != moving_pivot[0]) and (stopped_pivot[1] != moving_pivot[1])):
                    # Swap rows in the matrix
                    temp = matrix[stopped_pivot[0], :].copy()
                    matrix[stopped_pivot[0], :] = matrix[moving_pivot[0], :]
                    matrix[moving_pivot[0], :] = temp

                    # Update new_pivot_indices with the new positions
                    new_pivot_indices.remove(stopped_pivot)
                    new_pivot_indices.remove(moving_pivot)
                    new_pivot_indices.append((moving_pivot[0], stopped_pivot[1]))
                    new_pivot_indices.append((stopped_pivot[0], moving_pivot[1]))

                    swapped = True
                    break
            if swapped:
                break

        if not swapped:
            break

# Function to perform type I row operations to put pivots in order for inverse
def typeI_pivots_Inverse(matrix, identity_matrix, pivot_index):
    pivot_indices = pivot_index
    new_pivot_indices = pivot_indices
    while True:
        identity_row_indices, identity_column_indices = np.where(identity_matrix == 1)
        identity_pivot_indices = list(zip(identity_row_indices, identity_column_indices))

        swapped = False
        pivot_indices = new_pivot_indices

        for stopped_pivot in pivot_indices:
            for moving_pivot in pivot_indices:
                if ((stopped_pivot[0] > moving_pivot[0]) and (stopped_pivot[1] < moving_pivot[1])) or \
                ((stopped_pivot[0] < moving_pivot[0]) and (stopped_pivot[1] > moving_pivot[1])) and ((stopped_pivot[0] != moving_pivot[0]) and (stopped_pivot[1] != moving_pivot[1])):
                    # Swap rows in the matrix
                    temp = matrix[stopped_pivot[0], :].copy()
                    matrix[stopped_pivot[0], :] = matrix[moving_pivot[0], :]
                    matrix[moving_pivot[0], :] = temp

                    # Update new_pivot_indices with the new positions
                    new_pivot_indices.remove(stopped_pivot)
                    new_pivot_indices.remove(moving_pivot)
                    new_pivot_indices.append((moving_pivot[0], stopped_pivot[1]))
                    new_pivot_indices.append((stopped_pivot[0], moving_pivot[1]))

                    # Swap rows in the identity matrix
                    temp_identity = identity_matrix[stopped_pivot[0], :].copy()
                    identity_matrix[stopped_pivot[0], :] = identity_matrix[moving_pivot[0], :]
                    identity_matrix[moving_pivot[0], :] = temp_identity

                    swapped = True
                    break
            if swapped:
                break

        if not swapped:
            break

# Function to perform type I row operations to move zero rows to the bottom
def typeI_zeroRows_RREF(matrix):
    while True:
        rowOfZeros = np.all(matrix == 0, axis=1)
        ZeroRow_indices = np.where(rowOfZeros)[0]
        rowOfNonZeros = np.any(matrix != 0, axis=1)
        NonZeroRow_indices = np.where(rowOfNonZeros)[0]

        swapped = False
        for zero_index in ZeroRow_indices:
            for nonZero_index in NonZeroRow_indices:
                if ((zero_index < nonZero_index)):
                    # Perform the swap
                    temp = matrix[zero_index, :].copy()  
                    matrix[zero_index, :] = matrix[nonZero_index, :]
                    matrix[nonZero_index, :] = temp
                    swapped = True
                    break
            if swapped:
                break

        if not swapped:
            break

# Function to perform type I row operations to move zero rows to the bottom for inverse
def typeI_zeroRows_Inverse(matrix, identity_matrix):
    while True:
        rowOfZeros = np.all(matrix == 0, axis=1)
        ZeroRow_indices = np.where(rowOfZeros)[0]
        rowOfNonZeros = np.any(matrix != 0, axis=1)
        NonZeroRow_indices = np.where(rowOfNonZeros)[0]

        swapped = False
        for zero_index in ZeroRow_indices:
            for nonZero_index in NonZeroRow_indices:
                if ((zero_index < nonZero_index)):
                    # Perform the swap
                    temp = matrix[zero_index, :].copy()  
                    matrix[zero_index, :] = matrix[nonZero_index, :]
                    matrix[nonZero_index, :] = temp

                    temp_identity = identity_matrix[zero_index, :].copy()  
                    identity_matrix[zero_index, :] = identity_matrix[nonZero_index, :]
                    identity_matrix[nonZero_index, :] = temp_identity
                    swapped = True
                    break
            if swapped:
                break

        if not swapped:
            break

# Function to get the Reduced Row Echelon Form (RREF) of a matrix
def getRREF(matrix):
    copy_matrix = matrix.copy()

    pivot_index = typeII_and_III_RREF(matrix)
   
    # Type I: Put pivots in order
    typeI_pivots_RREF(matrix, pivot_index)
    # Type I: Move zero rows to the bottom
    typeI_zeroRows_RREF(matrix)

    print("Original Matrix: ")
    np.set_printoptions(precision=5, suppress=True)
    print(copy_matrix)

    # Printing the matrix using NumPy's print options
    print("Reduced Row Echelon Form (RREF):")
    np.set_printoptions(precision=5, suppress=True)
    matrix_string = matrix_to_string(matrix)
    print(matrix_string)
    
    return matrix

# Function to get the inverse of a matrix
def getInverse(matrix):
    numberOfRows, numberOfColumns = matrix.shape
    matrix_copy = matrix.copy()
    identity_matrix = np.identity((numberOfRows))

    if numberOfRows != numberOfColumns:
        print("Matrix must be square to compute inverse")
        quit()

    if(getDeterminant(matrix) == 0):
        print("Matrix is singular")
        quit()
    
    pivot_index = typeII_and_III_Inverse(matrix, identity_matrix)

    # Type I: Put pivots in order
    typeI_pivots_Inverse(matrix, identity_matrix, pivot_index)

    # Type I: Move zero rows to the bottom
    typeI_zeroRows_Inverse(matrix, identity_matrix)
    
    for (i, j), value in np.ndenumerate(matrix):
        matrix[i,j] = fix_near_zero(matrix[i,j])

    print("Original Matrix: ")
    np.set_printoptions(precision=10, suppress=True)
    matrix_copy_str = matrix_to_string(matrix_copy)
    print(matrix_copy_str)

    print("Inverse Matrix: ")
    np.set_printoptions(precision=10, suppress=True)
    identity_matrix_string = matrix_to_string(identity_matrix)
    print(identity_matrix_string)

    print("Reduced Row Echelon Form (RREF):")
    np.set_printoptions(precision=5, suppress=True)
    matrix_string = matrix_to_string(matrix)
    print(matrix_string)

# Function to get the minor of a matrix
def get_ijMinor(matrix, row, column):
    return np.delete(np.delete(matrix, row, axis=0), column, axis=1)

# Function to get the determinant of a matrix
def getDeterminant(matrix):
    numberOfRows, numberOfColumns = matrix.shape

    if numberOfRows != numberOfColumns:
        print("Matrix must be square to compute determinant")
        quit()

    if numberOfRows == 1:
        return matrix[0, 0]

    if numberOfRows == 2:
        return (matrix[0, 0] * matrix[1, 1]) - (matrix[0, 1] * matrix[1, 0])

    determinant = 0
    for col in range(numberOfColumns):
        sign = (-1) ** col
        minor = get_ijMinor(matrix, 0, col)
        cofactor = matrix[0, col]
        determinant += sign * cofactor * getDeterminant(minor)
    return determinant

# Function to start the program
def start():
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows * numberOfColumns

    entries_str = input("Enter all entries from right to left with a comma in between each one").split(',')

    if(len(entries_str) != size):
        print("Invalid number of entries")
        quit()
    else:
        entries_float = [float(entry) for entry in entries_str]
        matrix = np.array(entries_float).reshape(numberOfRows, numberOfColumns)

    getRREF(matrix)

# Function to get matrix from file
def fromFile(numberOfRows):
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows * numberOfColumns

    entries_str = input("Enter all entries from right to left with a comma in between each one").split(',')

    if(len(entries_str) != size):
        print("Invalid number of entries")
        quit()
    else:
        entries_float = [float(entry) for entry in entries_str]
        matrix = np.array(entries_float).reshape(numberOfRows, numberOfColumns)

    getRREF(matrix, numberOfRows)

# Function to generate a random matrix and get its RREF
def random():
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows * numberOfColumns

    positive_infinity = 1e9
    negative_infinity = -1e9  

    random_matrix = np.random.uniform(low=negative_infinity, high=positive_infinity, size=(numberOfRows, numberOfColumns))
    print(random_matrix)
    getRREF(random_matrix)

# Function to generate a random matrix and get its inverse
def random_inverse():
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows * numberOfColumns

    positive_infinity = 1e9
    negative_infinity = -1e9  

    random_matrix = np.random.uniform(low=negative_infinity, high=positive_infinity, size=(numberOfRows, numberOfColumns))
    print(random_matrix)
    getInverse(random_matrix)

# Function to generate a random matrix and get its determinant
def random_determinant():
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows * numberOfColumns

    positive_infinity = 1e9
    negative_infinity = -1e9  

    random_matrix = np.random.uniform(low=negative_infinity, high=positive_infinity, size=(numberOfRows, numberOfColumns))
    print(random_matrix)
    determinant = getDeterminant(random_matrix)
    print(f"Determinant is: {determinant}")

# Function to get the inverse of a matrix from user input
def getInverseStart():
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows * numberOfColumns

    positive_infinity = 1e6  
    negative_infinity = -1e6 

    entries_str = input("Enter all entries from right to left with a comma in between each one").split(',')

    if(len(entries_str) != size):
        print("Invalid number of entries")
        quit()
    else:
        entries_float = [float(entry) for entry in entries_str]
        matrix = np.array(entries_float).reshape(numberOfRows, numberOfColumns) 
        
        np.set_printoptions(precision=5, suppress=True)
        print(matrix)

    getInverse(matrix)


# Function to get the determinant of a matrix from user input
def getDeterminantStart():
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))

    if numberOfRows != numberOfColumns:
        print("Matrix must be square to compute determinant")
        return

    entries_str = input("Enter all entries from right to left with a comma in between each one (no spaces): ").split(',')

    if len(entries_str) != numberOfRows * numberOfColumns:
        print("Invalid number of entries")
        return

    entries_float = [float(entry) for entry in entries_str]
    matrix = np.array(entries_float).reshape(numberOfRows, numberOfColumns)
    
    determinant = getDeterminant(matrix)
    print(f"Determinant is: {determinant}")

# Function to run the program
def run_Program():
    userChoice_RREF = int(input("Choose which do you want to run: "
                       "\n [0] Get the RREF of a Matrix"
                       "\n [1] Get Inverse of a Matrix"
                       "\n [2] Get the Determinant of a Matrix"
                       "\n [3] Get the RREF of a Random Matrix"
                       "\n [4] Get the Inverse of a Random Matrix"
                       "\n [5] Get the Determinant of a Random Matrix\n"
                       ))

    numberOfRows = 0;

    if userChoice_RREF == 0:
        start()
    elif userChoice_RREF == 1:
        getInverseStart()
    elif userChoice_RREF == 2:
        getDeterminantStart()
    elif userChoice_RREF == 3:
        random()
    elif userChoice_RREF == 4:
        random_inverse()
    elif userChoice_RREF == 5:
        random_determinant()



run_Program()
    


