import numpy as np
import sympy as sp


def getRREF(matrix, numberOfRows):
    x = sp.symbols('x')
    pivots = []
    numberOfRows, numberOfColumns = matrix.shape
    
    # The first non-zero entry of each row must be 1
    for i in range(numberOfRows):
        first_nonZero = False  # Used to flag if we found the pivot
        for jth_column, entry in enumerate(matrix[i, :]):  # Traversing for every entry in the ith row (traversing columns of the ith row)
            if not first_nonZero:
                # If this is a non-zero element, divide the entire row by it to make it a leading 1 and mark that we found the pivot
                if entry != 0:
                    # Use numpy.where to find the indices of the entry
                    pivots.append((i, jth_column))
                    matrix[i, :] /= entry
                    first_nonZero = True

                    

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

    # Replace -0.00 with 0.00
    matrix[matrix == -0.0] = 0.0
    for stopped_pivot in pivots:
        for moving_pivot in pivots:
            if((stopped_pivot[0] > moving_pivot[0]) and (stopped_pivot[1] < moving_pivot[1])):
                temp = matrix[stopped_pivot[0], :].copy()
                matrix[stopped_pivot[0], :] = matrix[moving_pivot[0], :]
                matrix[moving_pivot[0], :]  = temp
        #I have a pivot with a 1 which has a greater row but a less than column than pivots ABOVE it swap)
    # Type I: Move zero rows to the bottom
    while True:
        rowOfZeros = np.all(matrix == 0, axis=1)
        ZeroRow_indices = np.where(rowOfZeros)[0]
        rowOfNonZeros = np.any(matrix != 0, axis=1)
        NonZeroRow_indices = np.where(rowOfNonZeros)[0]

        swapped = False
        for zero_index in ZeroRow_indices:
            for nonZero_index in NonZeroRow_indices:
                if zero_index < nonZero_index:
                    # Perform the swap
                    temp = matrix[zero_index, :].copy()  # Use copy to avoid referencing the same memory
                    matrix[zero_index, :] = matrix[nonZero_index, :]
                    matrix[nonZero_index, :] = temp
                    swapped = True
                    break  # Exit the inner loop to update indices
            if swapped:
                break  # Exit the outer loop to update indices

        if not swapped:
            break  # If no swaps were made, exit the while loop

    # Printing the matrix using NumPy's print options
    print("Reduced Row Echelon Form (RREF):")
    np.set_printoptions(precision=2, suppress=True)
    print(matrix)
    
    return matrix

        


       


def start(numberOfRows):

    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows*numberOfColumns



    #want to make this also work by taking direclty from files
    entries_str = list(input("Enter all entries from right to left with a comma in between each one").replace(',', ''))


    if(len(entries_str) != size):
        print("Invalid number of entires")
        quit()
    else:
        #List comprehenssion 
        #Turn all elements in list into ints
        entries_float = [float(entry) for entry in entries_str]

        matrix = np.array(entries_float).reshape(numberOfRows, numberOfColumns)

    getRREF(matrix, numberOfRows)

    def fromFile(numberOfRows):

        numberOfRows = int(input("What is the number of rows of the matrix?"))
        numberOfColumns = int(input("What is the number of columns of the matrix?"))
        size = numberOfRows*numberOfColumns



        #want to make this also work by taking direclty from files
        entries_str = list(input("Enter all entries from right to left with a comma in between each one").replace(',', ''))


        if(len(entries_str) != size):
            print("Invalid number of entires")
            quit()
        else:
            #List comprehenssion 
            #Turn all elements in list into ints
            entries_float = [float(entry) for entry in entries_str]

            matrix = np.array(entries_float).reshape(numberOfRows, numberOfColumns)

        getRREF(matrix, numberOfRows)

def random(numberOfRows):
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows*numberOfColumns
    
    positive_infinity = 1e6  
    negative_infinity = -1e6  

    random_matrix = np.random.uniform(low= negative_infinity, high= positive_infinity, size=(numberOfRows, numberOfColumns))
    
    getRREF(random_matrix, numberOfRows)


def run_Program():
    userChoice_RREF = int(input("Choose which do you want to run: "
                       "\n [0] Random Matrix "
                       "\n [1] Custom User Entered Matrix \n"))

    numberOfRows = 0;

    if userChoice_RREF == 0:
        random(numberOfRows)
    elif userChoice_RREF == 1:
        start(numberOfRows)
    else:
        start(numberOfRows)


run_Program()
