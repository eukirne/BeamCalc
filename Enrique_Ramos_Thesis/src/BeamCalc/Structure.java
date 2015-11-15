/**
 * (./) Structure.java v0.1 05/09/2011
 * @author Enrique Ramos Melgar
 * http://www.esc-studio.com
 *
 * THIS LIBRARY IS RELEASED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE
 * http://creativecommons.org/licenses/by/3.0/
 * http://www.processing.org/
 */

package BeamCalc;

import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.text.DecimalFormat;
import java.util.ArrayList;
import orbit.library.Orbit;
import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PFont;
import processing.core.PGraphics;
import processing.core.PVector;

/**
 * The class Structure contains is the main class in the library. All other
 * classes refer to this class, or are called by this class.
 * 
 * @author Enrique Ramos
 */
public class Structure {

	PApplet p; // Parent Processing Application
	Toolbox toolbox; // Frame for sliders and buttons
	BeamControl beamCont; // Frame for element control
	PDisp_Data s; // Additional Processing Frames
	ElementData eData;
	PGraphics selBuffer; // Selection Buffer for 3d Picking
	PFont fontA; // Font
	DXFImporter dxf; // Dxf Importer
	Orbit orb; // Orbit Object
	DecimalFormat df = new DecimalFormat("#0.000"); // Decimal Format
	Gradient grad = new Gradient(); // Gradient for visualisations

	// Main Node and Element Containers
	public ArrayList<Node> nodes;
	public ArrayList<Element> elms;

	// Focused Node and Element Containers
	public ArrayList<Node> focusNodes;
	public ArrayList<Element> focusElms;

	public int iterations = 10; // Number of iterations per drawn frame
	private int numPoints; // Number of intermediate points per element
	private int desiredFPS = 15; // Desired Frames per Second
	private int selNode = 0; // Index of selected node

	// Button Pressed modes for Mouse Behaviours
	private boolean shiftMode;
	private boolean controlMode;
	private boolean altMode;
	private boolean moving = false;

	public int prevMouseX, prevMouseY;
	protected int mode;
	boolean rendering;
	public boolean displayAxial;
	public float minAxial;
	public float maxAxial;
	float maxMoment;
	public boolean displayMoment;
	float minMoment;
	int displayMode;
	float minTorsion;
	float maxTorsion;
	float minDisp;
	float maxDisp;
	private int selBeam;
	private boolean hoverNode;
	private boolean hoverElm;
	private int intMode;
	private float fps;
	private float time;
	private float prevTime;
	boolean analysing = true;
	protected boolean drawing = true;
	public int topology;
	public boolean changeTopology = false;
	boolean graph = false;
	boolean record = false;

	/**
	 * Main constructor. It creates an empty instance of the Structure class.
	 **/
	public Structure(PApplet p, int numPoints) {
		this.p = p;
		this.numPoints = numPoints;

		// Methods for Library Management
		p.registerDispose(this);
		p.registerMouseEvent(this);
		p.registerKeyEvent(this);

		// Create the font in data folder
		fontA = p.createFont("ArialMT", 200);
		p.textFont(fontA, 200);

		// Initialise ArrayLists
		this.nodes = new ArrayList<Node>();
		this.elms = new ArrayList<Element>();
		this.focusNodes = new ArrayList<Node>();
		this.focusElms = new ArrayList<Element>();

		// Initialise Additional Frames
		toolbox = new Toolbox(this);
		beamCont = new BeamControl(this);
		s = new PDisp_Data(this);
		eData = new ElementData(this);

		// Initialise Orbit
		selBuffer = p.createGraphics(p.width, p.height, PApplet.P3D);
		orb = new Orbit(p, 0);
		orb.setCSScale(10000);
		orb.drawCS = true;

		// DXF Importer
		dxf = new DXFImporter(p, this);

		// Bring the main frame into focus
		p.frame.toFront();

	}

	/**
	 * This constructor creates an instance of structure from predefined node
	 * and element arrayLists
	 **/
	Structure(PApplet p, ArrayList<Node> inodes, ArrayList<Element> ibeams,
			int numPoints) {

		this(p, numPoints);

		this.nodes = new ArrayList<Node>(inodes);
		this.elms = new ArrayList<Element>(ibeams);

	}

	/**
	 * This method calls the evaluation functions of nodes and elements
	 **/
	public void eval() {

		if (analysing) {

			// Analyse the equilibrium equations for every node and element
			intMode = 1;

			for (int j = 0; j < iterations; j++) {

				switch (intMode) {
				case 0:
					euler(0.00001f, 0.0000000001f);
					break;
				case 1:
					verlet(0.001f, 0.000001f);
					break;
				}
			}

			// Analyse internal stresses in the elements
			for (int i = 0; i < elms.size(); i++) {
				elms.get(i).updateInternal();
			}

		}
	}

	/**
	 * This method evaluates an Euler step on every node in the structure.
	 **/
	private void euler(float h1, float h2) {

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).eval();
		}

		for (int i = 0; i < nodes.size(); i++) {
			nodes.get(i).euler(h1, h2);
		}

	}

	/**
	 * This method evaluates a forward-backward Euler step on every node in the
	 * structure.
	 **/
	private void verlet(float h1, float h2) {

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).eval();
		}

		for (int i = 0; i < nodes.size(); i++) {
			nodes.get(i).verlet(h1, h2);
		}

	}

	/**
	 * This method draws the current state of the structure on the window taken
	 * as an input. The beams currently focused will display stress graphs
	 **/
	public void draw(PApplet window) {

		if (drawing) {
			if (graph) {
				drawGraph(50, 50, 700, 500);
			} else {

				// Recalculate the maximum iterations for the
				// desired Frames per Second
				time = p.millis() - prevTime;
				prevTime = p.millis();
				fps = 1000 / time;
				iterations = (int) (1 + iterations * fps / desiredFPS);

				// Select display mode
				switch (displayMode) {
				case 0:
					break;
				case 1:
					showAxial();
					showScale("Axial Stresses (KN)", minAxial, maxAxial, "KN",
							0.001f);
					break;
				case 2:
					showMoment();
					showScale("Moment (KNm)", minMoment, maxMoment, "KNm",
							0.001f);
					break;
				case 3:
					showTorsion();
					showScale("Torsion (KNm)", minTorsion, maxTorsion, "KNm",
							0.001f);
					break;
				case 4:
					showDisp();
					showScale("Displacement (mm)", minDisp, maxDisp, "mm", 1f);
					break;
				}

				if (record == true) {
					p.beginRaw(PConstants.DXF, "output.dxf"); // Start recording
																// to the file
				}

				orb.pushOrbit(window); // Start Orbiting

				// Draws all the elements
				for (int i = 0; i < elms.size(); i++) {
					elms.get(i).draw();
				}

				// Draws all the nodes
				for (int i = 0; i < nodes.size(); i++) {
					nodes.get(i).draw(window);
				}

				orb.popOrbit(window); // Finish Orbiting

				if (record == true) {
					p.endRaw();
					record = false; // Stop recording to the file
				}

			}
		}
	}

	/**
	 * Create a box showing the scale of the values passed as parameters
	 * 
	 * @param title
	 * @param minVal
	 * @param maxVal
	 * @param units
	 * @param scale
	 */
	private void showScale(String title, float minVal, float maxVal,
			String units, float scale) {

		int rectScaleX = 15;
		int rectScaleY = 10;
		float numBoxes = 10;

		grad.setBounds(0, numBoxes);

		p.fill(0);
		p.noStroke();
		// p.stroke(0);
		p.textFont(fontA, 12);
		p.textAlign(PConstants.RIGHT);
		p.pushMatrix();
		p.translate(p.width - 80, 20);
		p.text(title, 50 + rectScaleX, 0);
		p.translate(0, 15);
		p.text("Min: " + df.format(minVal * scale) + " " + units,
				2 * rectScaleX, 10);
		for (int i = 0; i < numBoxes + 1; i++) {
			p.fill(grad.getColor(i)[0], grad.getColor(i)[1],
					grad.getColor(i)[2]);
			p.rect(50, i * (rectScaleY + 2), rectScaleX, rectScaleY);
		}
		p.fill(50);
		p.text("Max: " + df.format(maxVal * scale) + " " + units,
				2 * rectScaleX, (numBoxes + 1) * (rectScaleY + 2));
		p.popMatrix();

	}

	/**
	 * This calls the addBeam method from predefined planes
	 **/
	public void addBeam(PVector pt1, PVector pt2, int aMode, int bMode,
			int plane) {
		PVector planeVec = new PVector();
		if (plane == 0)
			planeVec.y = 1;
		if (plane == 1)
			planeVec.x = 1;
		if (plane == 2)
			planeVec.z = 1;
		addBeam(pt1, pt2, aMode, bMode, planeVec);
	}

	/**
	 * This is the routine to add an element. The other nodes in the structure
	 * are checked to avoid node duplication The degrees of freedom are defined
	 * from the calling function
	 **/
	public void addBeam(PVector pt1, PVector pt2, int aMode, int bMode,
			PVector planeVec) {

		Node a = new Node(p, this, pt1);
		Node b = new Node(p, this, pt2);

		// The following function checks if the node is duplicate. That is, it's
		// closer than nodeThreshold
		float nodeThreshold = 1f;
		boolean aDup = false;
		boolean bDup = false;

		for (int i = 0; i < nodes.size(); i++) {
			if (a.pos.dist(nodes.get(i).pos) < nodeThreshold) {
				aDup = true;
				a = nodes.get(i);
			}
			if (b.pos.dist(nodes.get(i).pos) < nodeThreshold) {
				bDup = true;
				b = nodes.get(i);
			}
		}

		if (!aDup) {
			nodes.add(a);
			nodes.get(nodes.size() - 1).id = nodes.size();
		}
		if (!bDup) {
			nodes.add(b);
			nodes.get(nodes.size() - 1).id = nodes.size();
		}

		// This modifies the node's constraints
		if (aMode == 1)
			a.linked = true;
		if (bMode == 1)
			b.linked = true;

		if (aMode == 2)
			a.clamped = true;
		if (bMode == 2)
			b.clamped = true;

		if (aMode == 3) {
			a.clamped = true;
			a.linked = true;
		}
		if (bMode == 3) {
			b.clamped = true;
			b.linked = true;
		}

		// Finally, add the new element to the elements arrayList
		Element elmTemp = new Element(p, this, a, b, numPoints, planeVec);
		elms.add(elmTemp);
		elms.get(elms.size() - 1).id = elms.size();
	}

	/**
	 * This function runs through the connected nodes, checks if they are
	 * isolated and removes beams and nodes in consequence.
	 * 
	 * @param thisNode
	 */
	private void removeNode(Node thisNode) {

		for (int i = 0; i < thisNode.nodeBeams.size(); i++) {
			Node otherNode = thisNode.nodeBeams.get(i).getOtherNode(thisNode);
			otherNode.nodeBeams.remove(thisNode.nodeBeams.get(i));
			if (otherNode.isDisconnected())
				nodes.remove(otherNode);
			focusNodes.remove(otherNode);
			elms.remove(thisNode.nodeBeams.get(i));
			focusElms.remove(thisNode.nodeBeams.get(i));
		}
		nodes.remove(thisNode);
		focusNodes.remove(thisNode);
	}

	/**
	 * This function runs through the connected nodes, checks if they are
	 * isolated and removes beams and nodes in consequence.
	 * 
	 * @param thisElm
	 */
	private void removeBeam(Element thisElm) {

		thisElm.a.nodeBeams.remove(thisElm);
		if (thisElm.a.isDisconnected())
			removeNode(thisElm.a);

		thisElm.b.nodeBeams.remove(thisElm);
		if (thisElm.b.isDisconnected())
			removeNode(thisElm.b);

		elms.remove(thisElm);
		focusElms.remove(thisElm);
	}

	private void checkHover() {

		selBuffer.beginDraw();
		orb.pushOrbit(selBuffer);

		selBuffer.background(255);

		// Draw the elements on the selection buffer
		selBuffer.strokeWeight(20);
		for (int i = nodes.size(); i < elms.size() + nodes.size(); i++) {
			selBuffer.stroke((int) (i / (255 * 255)), (int) (i / 255), i % 255);
			selBuffer.line(elms.get(i - nodes.size()).a.visualPos.x,
					elms.get(i - nodes.size()).a.visualPos.y,
					elms.get(i - nodes.size()).a.visualPos.z,
					elms.get(i - nodes.size()).b.visualPos.x,
					elms.get(i - nodes.size()).b.visualPos.y,
					elms.get(i - nodes.size()).b.visualPos.z);
		}

		// Draw the nodes on the selection buffer
		selBuffer.strokeWeight(10);
		for (int i = 0; i < nodes.size(); i++) {
			selBuffer.stroke((int) (i / (255 * 255)), (int) (i / 255), i % 255);
			selBuffer.pushMatrix();
			selBuffer.translate(nodes.get(i).visualPos.x,
					nodes.get(i).visualPos.y, nodes.get(i).visualPos.z);
			selBuffer.box(300);
			selBuffer.popMatrix();
		}

		selBuffer.loadPixels();
		int c = selBuffer.pixels[p.mouseY * p.width + p.mouseX];

		// int alpha = (c >> 24) & 0xFF;
		int red = (c >> 16) & 0xFF;
		int green = (c >> 8) & 0xFF;
		int blue = c & 0xFF;

		int selected = red * 255 * 255 + green * 255 + blue;

		for (int i = 0; i < nodes.size(); i++) {
			nodes.get(i).hover = false;
		}

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).hover = false;
		}

		hoverNode = false;
		hoverElm = false;

		if (selected < nodes.size()) {
			hoverNode = true;
			selNode = selected;
			nodes.get(selNode).drawGizmo(selBuffer, 100, 5, 3);
			nodes.get(selNode).hover = true;
		} else {
			if (selected < nodes.size() + elms.size()) {
				hoverElm = true;
				selBeam = elms.indexOf(elms.get(selected - nodes.size()));
				elms.get(selBeam).hover = true;
			} else {
			}
		}

		orb.popOrbit(selBuffer);
		selBuffer.endDraw();
	}

	private void checkGizmo() {

		selBuffer.beginDraw();
		orb.pushOrbit(selBuffer);

		selBuffer.background(255);

		if (selNode < nodes.size()) {
			nodes.get(selNode).drawGizmo(selBuffer, 500, 3, 3);
		}
		selBuffer.loadPixels();
		int c = selBuffer.pixels[p.mouseY * p.width + p.mouseX];

		// int alpha = (c >> 24) & 0xFF;
		int red = (c >> 16) & 0xFF;
		int green = (c >> 8) & 0xFF;
		int blue = c & 0xFF;

		if (red == 255 && green == 0 && blue == 0)
			nodes.get(selNode).activeGizmo(0);
		if (red == 0 && green == 255 && blue == 0)
			nodes.get(selNode).activeGizmo(1);
		if (red == 0 && green == 0 && blue == 255)
			nodes.get(selNode).activeGizmo(2);

		orb.popOrbit(selBuffer);
		selBuffer.endDraw();

	}

	private void focusNodes() {

		focusNodes = new ArrayList<Node>();
		focusElms = new ArrayList<Element>();

		focusNodes.add(nodes.get(selNode));

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).selected = false;
			elms.get(i).selected = false;
		}
		for (int i = 0; i < nodes.size(); i++) {
			nodes.get(i).selected = false;
			nodes.get(i).hover = false;
		}
		nodes.get(selNode).selected = true;

		for (int i = 0; i < nodes.get(selNode).nodeBeams.size(); i++) {
			focusNodes.add(nodes.get(selNode).nodeBeams.get(i).getOtherNode(
					nodes.get(selNode)));
			focusElms.add(nodes.get(selNode).nodeBeams.get(i));

			nodes.get(selNode).nodeBeams.get(i)
					.getOtherNode(nodes.get(selNode)).selected = true;
			nodes.get(selNode).nodeBeams.get(i).selected = true;
		}

	}

	private void focusBeam(Element focusedElm) {

		focusNodes = new ArrayList<Node>();
		focusElms = new ArrayList<Element>();

		focusNodes.add(focusedElm.a);
		focusElms.add(focusedElm);

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).selected = false;
			elms.get(i).selected = false;
		}
		for (int i = 0; i < nodes.size(); i++) {
			nodes.get(i).selected = false;
			nodes.get(i).hover = false;
		}

		focusedElm.a.selected = true;
		focusedElm.b.selected = true;

		focusedElm.selected = true;

	}

	public void mouseEvent(MouseEvent event) {

		checkHover();

		switch (event.getID()) {
		case MouseEvent.MOUSE_PRESSED:
			prevMouseX = p.mouseX;
			prevMouseY = p.mouseY;

			checkGizmo();

			break;
		case MouseEvent.MOUSE_RELEASED:

			break;
		case MouseEvent.MOUSE_CLICKED:
			// do something for mouse clicked

			if (shiftMode) {
				if (p.mouseButton == PConstants.LEFT) {
					if (hoverNode)
						nodes.get(selNode).clamped = !nodes.get(selNode).clamped;
					// if (hoverElm)elms.get(selBeam).clamped =
					// !elms.get(selBeam).clamped;
				}
				if (p.mouseButton == PConstants.RIGHT) {
					if (hoverNode)
						removeNode(nodes.get(selNode));
					if (hoverElm)
						removeBeam(elms.get(selBeam));
				}
			} else if (controlMode) {

			} else {
				if (p.mouseButton == PConstants.LEFT) {
					focusNodes = new ArrayList<Node>();
					focusElms = new ArrayList<Element>();
					for (int i = 0; i < elms.size(); i++) {
						elms.get(i).selected = false;
						elms.get(i).selected = false;
					}
					for (int i = 0; i < nodes.size(); i++) {
						nodes.get(i).selected = false;
						nodes.get(i).hover = false;
					}
					if (hoverNode) {
						focusNodes();
					}
					if (hoverElm)
						focusBeam(elms.get(selBeam));
				} else {
					nodes.get(selNode).linked = !nodes.get(selNode).linked;
				}
			}
			break;
		case MouseEvent.MOUSE_DRAGGED:

			if (shiftMode) {

				if (altMode) {
					nodes.get(selNode).modify(
							(p.mouseX - prevMouseX + p.mouseY - prevMouseY), 1);
				} else {
					nodes.get(selNode).modify(
							(p.mouseX - prevMouseX + p.mouseY - prevMouseY), 0);
				}
			}

			break;
		case MouseEvent.MOUSE_MOVED:
			// umm... forgot
			break;
		}
	}

	public void keyEvent(KeyEvent e) {

		if (e.isShiftDown()) {
			shiftMode = true;
		} else {
			shiftMode = false;
		}
		if (e.isControlDown()) {
			controlMode = true;
		} else {
			controlMode = false;
		}
		if (e.isAltDown()) {
			altMode = true;
		} else {
			altMode = false;
		}

		if (p.key == 'm') {
			moving = !moving;
		}

		if (p.key == 'a') {
			focusAll();
		}

		if (p.key == 'r') {
			render(true);
		}

		if (p.key == 'x') {
			displayMode = 1;
		}
		if (p.key == 'z') {
			displayMode = 2;
		}
		if (p.key == 'v') {
			displayMode = 3;
		}
		if (p.key == 'c') {
			displayMode = 4;
		}
		if (p.key == 'n') {
			resetStructure();
			analysing = false;
		}
		if (p.key == 's') {
			analysing = !analysing;
			for (int i = 0; i < elms.size(); i++) {
				elms.get(i).evalWeight();
			}
			iterations = 1;
		}

		if (p.key == '=') {
			p.save(p.dataPath("001.png"));
		}
		if (p.key == '0') { // Press R to save the file
			graph = !graph;
		}
		if (p.key == 'i') { // Press R to save the file
			record = true;
		}
	}

	private void resetStructure() {
		for (int i = 0; i < nodes.size(); i++) {
			nodes.get(i).pos = nodes.get(i).initialPos.get();
			nodes.get(i).dispVec = new PVector();
			nodes.get(i).rot = new PVector();
		}

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).reset();
		}
	}

	public void render(boolean mode) {

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).render = mode;

		}
		rendering = mode;
	}

	void focusAll() {
		focusNodes = new ArrayList<Node>(nodes);
		focusElms = new ArrayList<Element>(elms);

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).selected = true;
		}
		for (int i = 0; i < nodes.size(); i++) {
			nodes.get(i).selected = true;
		}
	}

	public void showAxial() {

		float axial = 0;
		minAxial = 1000000000;
		maxAxial = 0;

		for (int i = 0; i < elms.size(); i++) {
			axial = elms.get(i).axial;
			if (axial > maxAxial)
				maxAxial = axial;
			if (axial < minAxial)
				minAxial = axial;
		}

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).setColorRange(1, minAxial, maxAxial);
		}
	}

	public void showMoment() {

		float moment = 0;
		minMoment = 1000000000;
		maxMoment = 0;

		for (int i = 0; i < elms.size(); i++) {
			moment = elms.get(i).getMaxMoment();
			if (moment > maxMoment)
				maxMoment = moment;
		}

		for (int i = 0; i < elms.size(); i++) {
			moment = elms.get(i).getMaxMoment();
			if (moment < minMoment)
				minMoment = moment;
		}

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).setColorRange(2, minMoment, maxMoment);
		}
	}

	private void showTorsion() {
		float torsion = 0;
		minTorsion = 1000000000;
		maxTorsion = 0;

		for (int i = 0; i < elms.size(); i++) {
			torsion = Math.abs(elms.get(i).torsion);
			if (torsion > maxTorsion)
				maxTorsion = torsion;
			if (torsion < minTorsion)
				minTorsion = torsion;
		}

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).setColorRange(3, minTorsion, maxTorsion);
		}

	}

	private void showDisp() {

		float disp = 0;
		minDisp = 1000000000;
		maxDisp = 0;

		for (int i = 0; i < elms.size(); i++) {
			disp = elms.get(i).getMaxNodeDisp();
			if (disp > maxDisp)
				maxDisp = disp;
			if (disp < minDisp)
				minDisp = disp;
		}

		for (int i = 0; i < elms.size(); i++) {
			elms.get(i).setColorRange(4, minDisp, maxDisp);
		}

	}

	public void changeTopology() {
		changeTopology = true;
		iterations = 1;
		topology++;
		if (topology > 5)
			topology = 0;
	}

	/**
	 * Resets the Structure, emptying the Node and Element ArrayLists
	 */
	public void reset() {

		// Initialise ArrayLists
		this.nodes = new ArrayList<Node>();
		this.elms = new ArrayList<Element>();
		this.focusNodes = new ArrayList<Node>();
		this.focusElms = new ArrayList<Element>();

	}

	/**
	 * Draws a graph of the time evolution of the specified value
	 * 
	 * @param graphX
	 * @param graphY
	 * @param graphW
	 * @param graphH
	 */
	public void drawGraph(int graphX, int graphY, float graphW, float graphH) {

		p.pushMatrix();
		p.translate(0, p.height);
		p.scale(1, -1, 0);
		p.textFont(fontA, 12);

		p.stroke(0);
		p.fill(0);
		p.strokeWeight(2);
		p.line(graphX, graphY, graphX + graphW, graphY);
		p.line(graphX, graphY, graphX, graphY + graphH);

		try {
			float yScale = (float) graphH / focusNodes.get(0).getYscaleGraph();
			float xScale = (float) graphW / focusNodes.get(0).getXscaleGraph();

			float dispStep = focusNodes.get(0).getYscaleGraph() / 10;
			p.stroke(100);
			p.strokeWeight(1);
			for (int i = 0; i < focusNodes.get(0).getYscaleGraph() / dispStep; i++) {

				p.line(graphX, graphY + i * dispStep * yScale, graphX + graphW,
						graphY + i * dispStep * yScale);
				p.pushMatrix();
				p.translate(graphX - 40, graphY + i * dispStep * yScale);
				p.scale(1, -1, 1);
				p.text((i * dispStep), 0, 0);
				p.popMatrix();

			}

			for (int i = 0; i < focusNodes.get(0).graph.size() - 1; i++) {
				p.stroke(200, 0, 0);
				p.strokeWeight(2);
				p.line(graphX + focusNodes.get(0).graphTimes.get(i) * xScale,
						graphY - focusNodes.get(0).graph.get(i).z * yScale,
						graphX + focusNodes.get(0).graphTimes.get(i + 1)
								* xScale,
						graphY - focusNodes.get(0).graph.get(i + 1).z * yScale);

				p.stroke(100);
				p.strokeWeight(1);
				if ((int) (focusNodes.get(0).graphTimes.get(i) / 1000) != (int) (focusNodes
						.get(0).graphTimes.get(i + 1) / 1000)) {
					p.line(graphX + focusNodes.get(0).graphTimes.get(i)
							* xScale, graphY, graphX
							+ focusNodes.get(0).graphTimes.get(i) * xScale,
							graphY + graphH);
					p.pushMatrix();
					p.translate(graphX + focusNodes.get(0).graphTimes.get(i)
							* xScale - 5, graphY - 15);
					p.scale(1, -1, 1);
					p.text(1 + Math
							.round(focusNodes.get(0).graphTimes.get(i) / 1000),
							0, 0);
					p.popMatrix();
				}
			}

		} catch (Exception exc) {
			p.textFont(fontA, 15);
			p.translate(0, p.height);
			p.scale(1, -1, 1);
			p.text("No Node Selected", p.width / 2, p.height / 2);
		}
		p.popMatrix();
	}

	public void dispose() {
		// anything in here will be called automatically when
		// the parent applet shuts down. for instance, this might
		// shut down a thread used by this library.
	}

}
