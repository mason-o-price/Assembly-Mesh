import assemblyMesh as am

if __name__ == "__main__":
    maxTriangles = 300
    startingTriangles = 3           
    
    # Specify the directory and file name for the interaction matrix
    interactionMatrixFile = r'C:\Users\mason\OneDrive\Desktop\code\python\research\interactionMatrices\T13_interactions.csv'

    # Specify the directory and file name for the binding angles matrix
    bindingAnglesFile = r'C:\Users\mason\OneDrive\Desktop\code\python\research\bindingAngles\T13_angles.csv'
    
    # In case the input file paths above don't work, create a placeholer triangle and let the user upload new input files
    try:
        myMesh = am.initializeMesh(interactionMatrixFile, bindingAnglesFile, maxTriangles) 
        myApp = am.App(myMesh, startingTriangles, maxTriangles, interactionMatrixFile, bindingAnglesFile) 
        am.showMesh(myMesh, myApp.ax, startingTriangles)
    except: 
        baseTriangle = am.Triangle(index=0,  species=0, interactionMatrix=[], color=am.createColorMap(1))
        myMesh = [baseTriangle]
        myApp = am.App(myMesh, startingTriangles=1, maxTriangles=maxTriangles, interactionMatrixFile='', bindingAnglesFile='')
        am.showMesh(myMesh, myApp.ax, startingTriangles)
