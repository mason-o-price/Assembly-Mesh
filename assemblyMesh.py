# Import libraries
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import csv
import warnings
import time
from tkinter import Tk
from tkinter.filedialog import asksaveasfilename
from tkinter import filedialog as fd
from matplotlib.widgets import Button, Slider, TextBox
from matplotlib.colors import LinearSegmentedColormap

class Triangle:
    def __init__(self, index: int, species: int, interactionMatrix: list, color = 'b'):
        self.coordinates = np.array([[-0.5,0,0], [0.5,0,0], [0,np.sqrt(3)/2,0]]) # Matrix of vertex coordinates. # Format: [X | Y | Z], i.e., col1 = X, col2 = Y, col3 = Z, row1 = vertex1, row2 = vertex2, row3 = vertex3
        self.indx = index                               # Index. 0 to N-1 for N triangles.
        self.species = species                          # Species. 0 to k-1 for k species.
        self.color = color                              # Color. Format: [r, g, b].
        self.friends = [None for _ in range(3)]                               # Set of complementary species. Format: [[species1, side], [species2, side], [species3, side]]
        for side in [0, 1, 2]:
            for rowIdx, row in enumerate(interactionMatrix):
                if row[3*self.species + side] == 1:
                    self.friends[side] = ([rowIdx//3, rowIdx%3])
        self.boundEdges = [False, False, False]         # List of which edges are bound. Format: [edge 1 True/False, edge 2 True/False, edge 3 True/False]
        self.isFullyBound = False                       # Boolean storing whether each side has been bound

    def rotate(self, side: int, angle: float):
        if side not in [0,1,2]:
            Warning(f"Triangle index {self.indx}, species {self.species} \n -> Edge number must be 0,1 or 2. <-")

        # Define the order of the vertices based on the side number.
        vtxToMove = (side + 1)%3
        v1 = (side + 2)%3
        v2 = side
        
        # Calculate the edge vector connecting the two vertices on the given side
        edgeVector = np.array(self.coordinates[v2,:] - self.coordinates[v1,:])
        rotAxis = edgeVector/np.linalg.norm(edgeVector) # Normalize the rotation axis

        # Implement Rodriguez' rotation formula
        axisOrigin = (self.coordinates[v1,:] + self.coordinates[v2,:])/2    # Define the origin of the rotation axis to be the midpoint of the opposing edge
        XYZ0 = self.coordinates[vtxToMove,:] - axisOrigin   # Offset the coordinates of the vertex so that the axis origin is at [0,0,0]
        XYZrot = XYZ0*np.cos(angle) + np.cross(rotAxis, XYZ0)*np.sin(angle) + rotAxis*(np.dot(rotAxis, XYZ0))*(1-np.cos(angle)) # Apply the rotation
        XYZnew = XYZrot + axisOrigin                                        # Undo the initial offset

        # Update the coordinates of the vertex
        self.coordinates[vtxToMove,:] = XYZnew

    def updateBoundStatus(self):
        if self.boundEdges == [True, True, True]:
            self.isFullyBound = True

    def addNeighbor(self, side: int, mesh: list, numTriangles: int, interactionMatrix: list, bindingAngles: list, colorMatrix: list, vtxProximityThresh: float): 
        # Make sure the chosen side is a valid number
        if side not in [0, 1, 2]:
            Warning("Side number must be 0, 1, or 2.")

        # Make sure the base triangle is in the mesh to which we're adding neighbors.
        if not self in mesh:
            warnings.warn(f"Triangle {self.indx} (species {self.species}) not in mesh.")

        # If the triangle edge is already bound, skip it.
        if self.boundEdges[side]:
            return mesh, numTriangles

        # Find the complementary triangle species and side
        newSpecies, newTriangleSide = findComplementaryTriangle(self.species, side, interactionMatrix)

        # Create a new triangle
        newTriangle = Triangle(index = numTriangles, species = newSpecies, interactionMatrix = interactionMatrix, color = colorMatrix[newSpecies])

        # Update the number of triangles
        numTriangles += 1

        # Find the indices of the vertices that are/not shared in common between the two triangles.
        triangleACommonVtx1 = side                       # Index of the original triangle's 1st vertex shared in common with the new triangle. 0 -> 0, 1 -> 1, 2 -> 2
        triangleACommonVtx2 = (side + 2)%3               # Index of the original triangle's 2nd vertex shared in common with the new triangle. 0 -> 2, 1 -> 0, 2 -> 1
        triangleAUncommonVtx = (side + 1)%3              # Index of the original triangle's vertex that is not shared in common with the new triangle. 0 -> 1, 1 -> 2, 2 -> 0
        triangleBCommonVtx1 = (newTriangleSide - 1)%3    # Index of the new triangle's 1st vertex shared in common with the original triangle. 0 -> 2, 1 -> 0, 2 -> 1 
        triangleBCommonVtx2 = newTriangleSide            # Index of the new triangle's 2nd vertex shared in common with the original triangle. 0 -> 0, 1 -> 1, 2 -> 2 
        triangleBUncommonVtx = (newTriangleSide - 2)%3   # Index of the new triangle's vertex that is not shared in common with the original triangle. 0 -> 1, 1 -> 2, 2 -> 0

        # Determine the coordinates of the new triangle
        newTriangle.coordinates[triangleBCommonVtx1,:] = self.coordinates[triangleACommonVtx1,:]
        newTriangle.coordinates[triangleBCommonVtx2,:] = self.coordinates[triangleACommonVtx2,:]
        newTriangle.coordinates[triangleBUncommonVtx,:] = self.coordinates[triangleACommonVtx1,:] + self.coordinates[triangleACommonVtx2,:] - self.coordinates[triangleAUncommonVtx,:]

        # Rotate the new triangle about the bound edge by the binding angle
        newTriangle.rotate(side = newTriangleSide, angle = bindingAngles[self.species][newSpecies])

        # Update the status of the bound edges
        self.boundEdges[side] = True
        newTriangle.boundEdges[newTriangleSide] = True

        # Update the status of all edges being bound
        self.updateBoundStatus()
        newTriangle.updateBoundStatus()

        # Check if any of the edges overlap with the previous edges in the mesh (e.g., for self-closing assemblies).
        for (potentialSide, [friendlySpecies, friendlyEdge]) in enumerate(newTriangle.friends):
            potentialNeighbors = [mesh[i] for i,species in enumerate(findSpeciesVector(mesh)) if species == friendlySpecies]
            for triangleToTest in potentialNeighbors:
                # nearbyVertices = [dist for dist in np.linalg.norm(triangleToTest.coordinates - newTriangle.coordinates[triangleBUncommonVtx,:], axis = 1) if dist < vtxProximityThresh]
                if checkAdjacent(newTriangle, potentialSide, triangleToTest, friendlyEdge, vtxProximityThresh):  # If the triangles are bound
                    newTriangle.boundEdges[potentialSide] = True # Update the list of bound edges for triangle B
                    triangleToTest.boundEdges[friendlyEdge] = True # Update the list of bound edges for the complementary triangle

                    # Update the status of all edges being bound
                    newTriangle.updateBoundStatus()
                    triangleToTest.updateBoundStatus()

                    # Merge the vertices of the adjacent triangles
                    newTriangle, triangleToTest = mergeEdges(newTriangle, potentialSide, triangleToTest, friendlyEdge)

                    # Update the mesh
                    mesh[triangleToTest.indx] = triangleToTest

        # Update the mesh
        mesh[self.indx] = self
        mesh[newTriangle.indx] = newTriangle
        return mesh, numTriangles
    
def checkAdjacent(triangleA: Triangle, sideA: int, triangleB: Triangle, sideB: int, vtxProximityThresh: float) -> bool:
    # Find the vertex indices from the side numbers
    AVtx1 = sideA                                   # 0 -> 0, 1 -> 1, 2 -> 2
    AVtx2 = (sideA + 2)%3                           # 0 -> 2, 1 -> 0, 2 -> 1
    BVtx1 = (sideB - 1)%3                           # 0 -> 2, 1 -> 0, 2 -> 1 
    BVtx2 = sideB                                   # 0 -> 0, 1 -> 1, 2 -> 2

    dist1 =  np.linalg.norm(triangleA.coordinates[AVtx1,:] - triangleB.coordinates[BVtx1,:])
    dist2 =  np.linalg.norm(triangleA.coordinates[AVtx2,:] - triangleB.coordinates[BVtx2,:])

    if dist1<vtxProximityThresh and dist2<vtxProximityThresh:
        return True
    else:
        return False
    
def mergeEdges(triangleA: Triangle, sideA: int, triangleB: Triangle, sideB: int) -> Triangle:
    # Find the vertex indices from the side numbers
    AVtx1 = sideA                                   # 0 -> 0, 1 -> 1, 2 -> 2
    AVtx2 = (sideA + 2)%3                           # 0 -> 2, 1 -> 0, 2 -> 1
    BVtx1 = (sideB - 1)%3                           # 0 -> 2, 1 -> 0, 2 -> 1 
    BVtx2 = sideB                                   # 0 -> 0, 1 -> 1, 2 -> 2

    # Replace the coordinates of the vertices to be the average 
    triangleA.coordinates[AVtx1, :] = (triangleA.coordinates[AVtx1, :] + triangleB.coordinates[BVtx1, :])/2
    triangleB.coordinates[BVtx1, :] = (triangleA.coordinates[AVtx1, :] + triangleB.coordinates[BVtx1, :])/2
    triangleA.coordinates[AVtx2, :] = (triangleA.coordinates[AVtx2, :] + triangleB.coordinates[BVtx2, :])/2
    triangleB.coordinates[BVtx2, :] = (triangleA.coordinates[AVtx2, :] + triangleB.coordinates[BVtx2, :])/2

    return triangleA, triangleB


def findComplementaryTriangle(species: int, side: int, interactionMatrix: list) -> tuple:
    complementaryRow = [rowIdx for rowIdx, row in enumerate(interactionMatrix) if row[3*species+side]==1][0]
    complementarySpecies = complementaryRow//3
    complementarySide = complementaryRow%3
    return (complementarySpecies, complementarySide)

def findSpeciesVector(mesh: list[Triangle]) -> list:
    return [tri.species for tri in mesh if not tri is None]

def showMesh(mesh: list[Triangle], ax, trianglesToShow: int) -> None:
    # Retrieve the current view
    # dist = ax.get_proj()[3] if hasattr(ax, 'get_proj') else None
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    zlim = ax.get_zlim()


    # Set axis propertes
    ax.cla()                    # Clear axes
    ax.set_proj_type('ortho')   # Set orthogonal perspective
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_zlim(zlim)
    ax.set_aspect('equal')

    # Retrieve the coordinates and colors of the triangles
    XYZ = [tri.coordinates for tri in mesh[:trianglesToShow] if not tri is None]
    triColors = [tri.color for tri in mesh[:trianglesToShow] if not tri is None]

    ax.add_collection(Poly3DCollection(XYZ, facecolors = triColors, edgecolors='black'))

    # Format plot
    ax.grid(False)                                  # Hide grid lines
    ax.set_xticks([])                               # Hide x-axis ticks
    ax.set_yticks([])                               # Hide y-axis ticks
    ax.set_zticks([])                               # Hide z-axis ticks
    ax.set_axis_off()                               # Hide axis-lines
    
    plt.show()                                      # Show the plot

def generateMesh(baseMesh: list, numTriangles: int, interactionMatrix: list, bindingAngles: list, colorMatrix: list, maxTriangles: int) -> list:
    # Check that the size of the interaction matrix and binding angles matrix are compatible
    if len(interactionMatrix)/3 != len(bindingAngles):
        warnings.warn(f"The sizes of the interaction matrix and binding angle matrix do not match. \n Interaction Matrix: {len(interactionMatrix)} \n Binding angles: {len(bindingAngles)} \n -> Interaction matrix should be 3x the size of the binding angle matrix. <-")

    # Add triangles until we hit the max number of triangles (or we time out)
    timeStart = time.time()
    while time.time() - timeStart < 120:
        for triangleIdx in range(numTriangles):
            for side in [0, 1, 2]:
                if numTriangles >= maxTriangles:
                    print(f"{numTriangles} triangles")
                    return baseMesh
                if isClosed(baseMesh):
                    print(f"self-closure at {numTriangles} triangles")
                    return baseMesh
                tri = baseMesh[triangleIdx]
                if not tri.isFullyBound: # Check that we're not adding to a bound triangle
                    baseMesh, numTriangles = tri.addNeighbor(side, baseMesh, numTriangles, interactionMatrix, bindingAngles, colorMatrix, vtxProximityThresh=0.2)      # Add new triangles, and update existing ones.
                
    # If we've made it this far without returning something, the while loop timed out. 
    print(f"Timed out")
    return baseMesh

def isClosed(mesh: list) -> bool:
    closed = True
    for tri in mesh:
        if not tri is None:
            if not tri.isFullyBound:
                closed = False
    return closed

def createColorMap(numColors) -> list:
    # cmap = cm.roma
    videColors = np.array([[0,129,167], [0, 175, 185], [253, 252, 220], [254, 217, 183], [240, 113, 103]])/255
    cmap = LinearSegmentedColormap.from_list("myCmap", videColors, N=256)
    colors = cmap(np.linspace(0, 1, numColors))[:, :3]  # Ignore alpha channel if present
    return colors

def initializeMesh(interactionMatrixFile: str, bindingAnglesFile: str, maxTriangles: int) -> list:
    # Open the interaction matrix file
    with open(interactionMatrixFile, mode ='r') as intFile:
        reader = csv.reader(intFile)
        interactionMatrix = [[int(x) for x in line] for line in reader]

    # Open the binding angles matrix
    with open(bindingAnglesFile, mode ='r') as angFile:
        reader = csv.reader(angFile)
        bindingAngles = [[np.pi/180*float(x) for x in line] for line in reader]

    # Find the number of unique species
    numSpecies = len(bindingAngles)

    # Define a color matrix
    colorMatrix = createColorMap(numSpecies)

    # Create a triangle to seed the mesh
    baseTriangle = Triangle(index=0,
                            species=0,
                            interactionMatrix=interactionMatrix,
                            color=colorMatrix[0])

    # Create a mesh
    myMesh = [None]*maxTriangles
    myMesh[0] = baseTriangle                            # Initialize a mesh
    numTriangles = 1
    myMesh = generateMesh(myMesh, numTriangles, interactionMatrix, bindingAngles, colorMatrix, maxTriangles) # Add triangles to the mesh
    return myMesh

class App:
    def __init__(self, myMesh, startingTriangles, maxTriangles, interactionMatrixFile, bindingAnglesFile):
        self.myMesh = myMesh
        self.maxTriangles = maxTriangles
        # Create the figure
        self.fig = plt.figure("Assembly Mesh", figsize=(8, 4))
        plt.subplots_adjust(left=0.1, right=1, bottom=0, top=1)

        # Make a horizontal slider to control the number of triangles.
        self.axNum = self.fig.add_axes([0.135, 0.1, 0.35, 0.03])  # Left, bottom, width, height
        self.numSlider = Slider(
            ax=self.axNum,
            label='Number of  \n triangles  ',
            valmin=1,
            valmax=maxTriangles,
            valinit=startingTriangles,
            valstep=1,
        )

        # Create labels for the interaction matrix and binding angle files
        plt.figtext(0.025, 0.94, 'Interaction matrix')
        plt.figtext(0.025, 0.78, 'Binding angles')
        
        # Create axes for the text boxes storing the file locations
        axTextBoxInteraction = plt.axes([0.12, 0.86, 0.45, 0.06])  # Left, bottom, width, height
        axTextBoxAngles = plt.axes([0.12, 0.7, 0.45, 0.06])  # Left, bottom, width, height

        # Create text boxes for the files
        self.interactionTextBox = TextBox(axTextBoxInteraction, label = '', initial=interactionMatrixFile, color='white')
        self.anglesTextBox = TextBox(axTextBoxAngles, label = '', initial=bindingAnglesFile, color='white')

        # Create a `matplotlib.widgets.Button` to prompt the user to save the figure.
        buttonAxSave = self.fig.add_axes([0.8, 0.025, 0.125, 0.05])  # Left, bottom, width, height
        self.saveButton = Button(buttonAxSave, 'Save figure', hovercolor='0.975')
        self.saveButton.on_clicked(self.handleSaveButton)

        # Create a `matplotlib.widgets.Button` to prompt the user to import an interaction matrix
        buttonAxInt = self.fig.add_axes([0.025, 0.86, 0.075, 0.06])  # Left, bottom, width, height
        self.importButton1 = Button(buttonAxInt, 'Import', hovercolor='0.975')
        self.importButton1.on_clicked(self.handleImportInteractionButton)

        # Create a `matplotlib.widgets.Button` to prompt the user to import a binding angles matrix
        buttonAxAng = self.fig.add_axes([0.025, 0.7, 0.075, 0.06])  # Left, bottom, width, height
        self.importButton2 = Button(buttonAxAng, 'Import', hovercolor='0.975')
        self.importButton2.on_clicked(self.handleImportAnglesButton)

        # Create a `matplotlib.widgets.Button` to update the input information
        buttonAxUpdate = self.fig.add_axes([0.025, 0.6, 0.15, 0.06])  # Left, bottom, width, height
        self.buttonUpdate = Button(buttonAxUpdate, 'Update inputs', hovercolor='0.975')
        self.buttonUpdate.on_clicked(self.handleUpdateButton)
        
        # Create a text box to specify the max number of triangles
        axTextBoxMax = self.fig.add_axes([0.15, 0.25, 0.075, 0.06])  # Left, bottom, width, height
        # Create text boxes for the files
        self.MaxTextBox = TextBox(axTextBoxMax, label = 'Max triangles ', initial=str(maxTriangles))
        self.MaxTextBox.on_submit(self.handleMaxChange)

        # Prepare to plot 3D
        self.ax = self.fig.add_subplot(1, 2, 2, projection = '3d')
        self.ax.set_xlim(-3, 3)
        self.ax.set_ylim(-3, 3)
        self.ax.set_zlim(-1, 5)

        # Register the update function with the slider
        self.numSlider.on_changed(self.handleUpdateSlider)

    # The function to be called anytime a slider's value changes
    def handleUpdateSlider(self, val):
        setTriangles = self.numSlider.val
        showMesh(self.myMesh, self.ax, setTriangles)

    def handleSaveButton(self, event):
        # Hide the root window (we only want the dialog)
        root = Tk()
        root.withdraw()

        # Show save file dialog
        file_path = asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG (.png)", "*.png"), ("SVG (.svg)", "*.svg")],
            title="Save figure as..."
        )

        # If the user didn't cancel, save the figure
        if file_path:
            # Only save the bounding box of the target_ax
            extent = self.ax.get_window_extent().transformed(self.fig.dpi_scale_trans.inverted())
            self.fig.savefig(file_path, bbox_inches=extent)
            print(f"Figure saved to {file_path}")
        else:
            print("Save cancelled.")

    def handleImportInteractionButton(self, event):
        # Hide the root window (we only want the dialog)
        root = Tk()
        root.withdraw()
        
        filetypes = [('CSV files', '*.csv')]
        filename = fd.askopenfilename(
            title='Open file',
            initialdir='/',
            filetypes=filetypes)
        print(filename)

        self.interactionTextBox.set_val(filename)

    def handleImportAnglesButton(self, event):
        # Hide the root window (we only want the dialog)
        root = Tk()
        root.withdraw()
        
        filetypes = [('CSV files', '*.csv')]
        filename = fd.askopenfilename(
            title='Open file',
            initialdir='/',
            filetypes=filetypes)
        print(filename)

        self.anglesTextBox.set_val(filename)

    def handleUpdateButton(self, event):
        # Find the interaction matrix file
        interactionMatrixFile = self.interactionTextBox.text

        # Find the binding angles matrix
        bindingAnglesFile = self.anglesTextBox.text

        self.myMesh = initializeMesh(interactionMatrixFile, bindingAnglesFile, self.maxTriangles)

        # Show the mesh
        showMesh(self.myMesh, self.ax, self.numSlider.val)

        plt.show()

    def handleMaxChange(self, expression):
        # Find the specified new max number of triangles
        newMax = int(expression)

        # Check that we're actually changing something
        if newMax == self.maxTriangles:
            return
        # Update the properties of the app
        self.maxTriangles = newMax      

        # Find the current number of triangles to display
        currentNum = self.numSlider.val

        # Disable the previous slider
        self.numSlider.ax.remove()

        self.axNum = self.fig.add_axes([0.125, 0.1, 0.35, 0.03])  # Left, bottom, width, height

        # Update the slider for the number of triangles.
        self.numSlider = Slider(
            ax=self.axNum,
            label='Number of  \n triangles  ',
            valmin=1,
            valmax=newMax,
            valinit=currentNum,
            valstep=1,
        )

        # Find the interaction matrix file
        interactionMatrixFile = self.interactionTextBox.text

        # Find the binding angles matrix
        bindingAnglesFile = self.anglesTextBox.text

        # Initialize a new mesh
        self.myMesh = initializeMesh(interactionMatrixFile, bindingAnglesFile, self.maxTriangles)

        # Register the update function with the slider
        self.numSlider.on_changed(self.handleUpdateSlider)


        
