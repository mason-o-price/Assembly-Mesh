import matplotlib.pyplot as plt
import assemblyMesh as am

if __name__ == "__main__":
    # Specify the maximum number of triangles to generate for the mesh
    maxTriangles = 300

    # Specify the number of triangles to start with
    startingTriangles = 3
    
    # Specify the directory and file name for the interaction matrix
    interactionMatrixFile = r'C:\Users\mason\OneDrive\Desktop\code\python\research\interactionMatrices\T13_interactions.csv'

    # Specify the directory and file name for the binding angles matrix
    bindingAnglesFile = r'C:\Users\mason\OneDrive\Desktop\code\python\research\bindingAngles\T13_angles.csv'
    
    # Try to start the app with the input files from above
    try:
        myMesh = am.initializeMesh(interactionMatrixFile, bindingAnglesFile, maxTriangles) # Create the mesh
        myApp = am.App(myMesh, startingTriangles, maxTriangles, interactionMatrixFile, bindingAnglesFile) # Create the user interface
        am.showMesh(myMesh, myApp.ax, startingTriangles)                                   # Show the mesh
        plt.show()
    except Exception as ex: # In case the input file paths above don't work, create a placeholer triangle and let the user upload new input files
        print(f"Exception: {ex}")
        # Create a triangle to seed the mesh
        baseTriangle = am.Triangle(index=0,
                            species=0,
                            interactionMatrix=[],
                            color=am.createColorMap(1))
        myMesh = [baseTriangle]
        myApp = am.App(myMesh, 1, maxTriangles, '', '')
        am.showMesh(myMesh, myApp.ax, startingTriangles)
        plt.show()