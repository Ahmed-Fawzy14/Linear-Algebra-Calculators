import numpy as np
import sympy as sp



def getRREF(matrix):

    x = sp.symbols('x')
    #The first non-zero entry of each row must be 1
    for i in range(0,matrix.shape[0], 1):

        first_nonZero = False #Used to flag if we found the pivot

        for jth_column, entry in enumerate(matrix[i,:] ): #traversing for every entry in the ith row (traversing columns of the ith row)
            
            #If we have not reached the leading one continue
            if(first_nonZero == False): 
                '''If this is a non-zero element divide the entire row by it
                to make it a leading 1 and mark that we foudn the pivot'''
                if(entry != 0):
                    matrix[i,:]  /= entry
                    first_nonZero = True
                    for ith_row,column_entry in enumerate(matrix[:,jth_column]): #tranvesing in the column to kill
                        if((column_entry.size > 0) and (column_entry != 0) and (column_entry != matrix[i][jth_column])):
                            #change this to doing it using the actual ero so we can store it and figure out what happened to the determinant
                            #btw this is wrong you need to apply it to the whole row
                            #from this we can store the type 2 we used before this
                            #and type 3 inside this
                            #idt we need to use type 1
                            #if we have the types we can get the determinant
                            #if we have the types we can apply them to the identity to get the inverse if the det != 0
                            equation = sp.Eq(entry + x*matrix[ith_row][jth_column], 0)
                            solution = sp.solve(equation, x)
                            float_solution = [float(sol) for sol in solution]
                            matrix[:,jth_column] = matrix[:,jth_column] + float_solution*matrix[ith_row][jth_column]

                else:
                    pass
            else:
                break


    # Replace -0.00 with 0.00
    matrix[matrix == -0.0] = 0.0

  # Printing the matrix using NumPy's print options
    print("Reduced Row Echelon Form (RREF):")
    np.set_printoptions(precision=2, suppress=True)
    print(matrix)

def start():

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

    getRREF(matrix)

def random():
    numberOfRows = int(input("What is the number of rows of the matrix?"))
    numberOfColumns = int(input("What is the number of columns of the matrix?"))
    size = numberOfRows*numberOfColumns
    
    positive_infinity = 1e6  
    negative_infinity = -1e6  

    random_matrix = np.random.uniform(low= negative_infinity, high= positive_infinity, size=(numberOfRows, numberOfColumns))
    
    getRREF(random_matrix)


def run_Program():
    userChoice_RREF = int(input("Choose which do you want to run: "
                       "\n [0] Random Matrix "
                       "\n [1] Custom User Entered Matrix \n"))

    if userChoice_RREF == 0:
        random()
    elif userChoice_RREF == 1:
        start()
    else:
        start()


run_Program()
