import assemblyMesh as am

def main():
    maxTriangles = 300
    startingTriangles = 3           
    
    # Specify the directory and file name for the interaction matrix
    interactionMatrixFile = r''

    # Specify the directory and file name for the binding angles matrix
    bindingAnglesFile = r''
    
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

if __name__ == "__main__":
    main()
