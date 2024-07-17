import numpy as np
import sympy as sp


def matrix_to_string(matrix):
    matrix_str = "["
    for i, row in enumerate(matrix):
        row_str = " ".join(f"{int(x)}" if x.is_integer() else f"{x}" for x in row)
        matrix_str += f"[{row_str}]"
        if i < matrix.shape[0] - 1:
            matrix_str += "\n "
    matrix_str += "]"
    return matrix_str

def fix_near_zero(value, tolerance=1e-4):
    if np.isclose(value, 0 , atol=tolerance):
        return 0.0
    elif np.isclose(value, 1, atol=tolerance):
        return 1.0
    return value

def typeII_and_III_RREF(matrix):
    x = sp.symbols('x')
    pivot_indices = []
    numberOfRows, numberOfColumns = matrix.shape

    # The first non-zero entry of each row must be 1
    for i in range(numberOfRows):
        first_nonZero = False  # Used to flag if we found the pivot
        for jth_column, entry in enumerate(matrix[i, :]):  # Traversing for every entry in the ith row (traversing columns of the ith row)
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
                    for ith_row, column_entry in enumerate(matrix[:, jth_column]):  # Traversing in the column to eliminate other entries
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


def typeII_and_III_Inverse(matrix,identity_matrix):
    x = sp.symbols('x')
    pivot_indices = []
    numberOfRows, numberOfColumns = matrix.shape
    # The first non-zero entry of each row must be 1
    for i in range(numberOfRows):
        first_nonZero = False  # Used to flag if we found the pivot
        for jth_column, entry in enumerate(matrix[i, :]):  # Traversing for every entry in the ith row (traversing columns of the ith row)
            if not first_nonZero:
                
                # If this is a non-zero element, divide the entire row by it to make it a leading 1 and mark that we found the pivot
                #Fix near zero errors
                entry = fix_near_zero(entry)
                if entry != 0:
                    # Use numpy.where to find the indices of the entry
                    print(f"entry: {entry}")
                    print(f"matrix[i,:]: {matrix[i, :]}")
                    matrix[i, :] /= entry
                    identity_matrix[i, :] /= entry
                    first_nonZero = True
                    pivot_row_index = i 
                    pivot_column_index = jth_column
                    pivot_indices.append((pivot_row_index, pivot_column_index))

                    

                    # Type III: Make the other entries in this column zero
                    for ith_row, column_entry in enumerate(matrix[:, jth_column]):  # Traversing in the column to eliminate other entries
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


def typeI_pivots_RREF(matrix, pivot_index):
    pivot_indices = pivot_index  # Initialize pivot_indices to pivot_index initially
    new_pivot_indices = pivot_indices
    while True:
       

        swapped = False

        pivot_indices = new_pivot_indices  # Update pivot_indices at the end of each iteration

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
                break  # Exit the outer loop to update indices

        if not swapped:
            break  # If no swaps were made, exit the while loop

def typeI_pivots_Inverse(matrix, identity_matrix, pivot_index):
    pivot_indices = pivot_index  # Initialize pivot_indices to pivot_index initially
    new_pivot_indices = pivot_indices
    while True:
        identity_row_indices, identity_column_indices = np.where(identity_matrix == 1)
        identity_pivot_indices = list(zip(identity_row_indices, identity_column_indices))

        swapped = False

        pivot_indices = new_pivot_indices  # Update pivot_indices at the end of each iteration

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
                break  # Exit the outer loop to update indices

        if not swapped:
            break  # If no swaps were made, exit the while loop



        
def typeI_zeroRows_RREF(matrix):
     while True:
        #I need to add something like this in the top loop to update things
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
                    break  # Exit the inner loop to update indices
            if swapped:
                break  # Exit the outer loop to update indices

        if not swapped:
            break  # If no swaps were made, exit the while loop

def typeI_zeroRows_Inverse(matrix, identity_matrix):
    # Type I: Move zero rows to the bottom
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
                    break  # Exit the inner loop to update indices
            if swapped:
                break  # Exit the outer loop to update indices

        if not swapped:
            break  # If no swaps were made, exit the while loop



def getRREF(matrix):
    copy_matrix = matrix.copy()


    pivot_index = typeII_and_III_RREF(matrix)
   
    #Type I: Put pivots in order
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

   
   #Type I: Put pivots in order
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
    
def get_ijMinor(matrix, row, column):
    return np.delete(np.delete(matrix, row, axis=0), column, axis=1)

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




'''def getDeterminant(matrix, cols, minor_copy, counter):
    numberOfRows, numberOfColumns = matrix.shape

    if (counter >= 1):
        minor = get_ijMinor(minor_copy, 0, cols)
    else:
        minor = matrix

    minor_copy = minor
 
    numberOfRows_minor, numberOfColumns_minor = minor.shape
    if((numberOfRows_minor == 1) and (numberOfColumns_minor == 1)):
        return minor[0,0]
    
    if((numberOfRows_minor == 2) and (numberOfColumns_minor == 2)):
        det2x2 = (minor[0,0]*minor[1,1]) - (minor[0,1]*minor[1,0])
        return det2x2 #Determinant of 2x2
    
    else:
        #Sends matrix after minor
        counter = counter + 1
        det = getDeterminant(matrix, cols, minor_copy, counter)
        print(f"multiplying {minor[0,cols]} by {det} ")
        return minor[0,cols]*det'''
    

def start():

    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows*numberOfColumns



    #want to make this also work by taking direclty from files
    entries_str = input("Enter all entries from right to left with a comma in between each one").split(',')
    print(entries_str)


    if(len(entries_str) != size):
        print("Invalid number of entires")
        quit()
    else:
        #List comprehenssion 
        #Turn all elements in list into ints
        entries_float = [float(entry) for entry in entries_str]

        matrix = np.array(entries_float).reshape(numberOfRows, numberOfColumns)

    getRREF(matrix)


    def fromFile(numberOfRows):

        numberOfRows = int(input("What is the number of rows of the matrix?"))
        numberOfColumns = int(input("What is the number of columns of the matrix?"))
        size = numberOfRows*numberOfColumns



        #want to make this also work by taking direclty from files
        entries_str = input("Enter all entries from right to left with a comma in between each one").split(',')
        print(entries_str)
        #negatives don't work idk why?
        if(len(entries_str) != size):
            print("Invalid number of entires")
            quit()
        else:
            #List comprehenssion 
            #Turn all elements in list into ints
            entries_float = [float(entry) for entry in entries_str]

            matrix = np.array(entries_float).reshape(numberOfRows, numberOfColumns)

        getRREF(matrix, numberOfRows)

def random():
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows*numberOfColumns
    
    #At 1e4 or 1e3 it is good else we start to deal with large value errors 
    positive_infinity = 1e9
    negative_infinity = -1e9  

    random_matrix = np.random.uniform(low= negative_infinity, high= positive_infinity, size=(numberOfRows, numberOfColumns))
    print(random_matrix)
    getRREF(random_matrix)

def random_inverse():
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows*numberOfColumns

    #At 1e4 or 1e3 it is good else we start to deal with large value errors 
    positive_infinity = 1e9
    negative_infinity = -1e9  

    random_matrix = np.random.uniform(low= negative_infinity, high= positive_infinity, size=(numberOfRows, numberOfColumns))
    print(random_matrix)
    getInverse(random_matrix)

def random_determinant():
        numberOfRows = int(input("What is the number of rows of the matrix?"))
        numberOfColumns = int(input("What is the number of columns of the matrix?"))
        size = numberOfRows*numberOfColumns

        #At 1e4 or 1e3 it is good else we start to deal with large value errors 
        positive_infinity = 1e9
        negative_infinity = -1e9  

        random_matrix = np.random.uniform(low= negative_infinity, high= positive_infinity, size=(numberOfRows, numberOfColumns))
        print(random_matrix)
        determinant = getDeterminant(random_matrix)
        print(f"Determinant is: {determinant}")



def getInverseStart():
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows*numberOfColumns


    positive_infinity = 1e6  
    negative_infinity = -1e6 

    entries_str = input("Enter all entries from right to left with a comma in between each one").split(',')
    print(entries_str)
    #negatives don't work idk why?
    if(len(entries_str) != size):
        print("Invalid number of entires")
        quit()
    else:
        #RREF and so Inverse works for singular digits but not decimals and more than singular digits
        entries_float = [float(entry) for entry in entries_str]

        matrix = np.array(entries_float).reshape(numberOfRows, numberOfColumns) 
        
        np.set_printoptions(precision=5, suppress=True)
        print(matrix)

    getInverse(matrix)

def getDeterminantStart():
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))

    if numberOfRows != numberOfColumns:
        print("Matrix must be square to compute determinant")
        return

    entries_str = input("Enter all entries from right to left with a comma in between each one: ").split(',')
    if len(entries_str) != numberOfRows * numberOfColumns:
        print("Invalid number of entries")
        return

    entries_float = [float(entry) for entry in entries_str]
    matrix = np.array(entries_float).reshape(numberOfRows, numberOfColumns)
    
    determinant = getDeterminant(matrix)
    print(f"Determinant is: {determinant}")

def run_Program():
    userChoice_RREF = int(input("Choose which do you want to run: "
                       "\n [0] Get the RREF of a Matrix"
                       "\n [1] Get Inverse of a Matrix"
                       "\n [2] Get the Determinant of a Matrix"
                       "\n [3] Get the RREF of a Random Matrix"
                       "\n [4] Get the Inverse of a Random Matrix"
                       "\n [5] Get the Determinant of a Random Matrix \n"
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