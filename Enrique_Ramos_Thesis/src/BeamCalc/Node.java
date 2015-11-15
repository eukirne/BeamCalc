/**
 * (./) Node.java v0.1 05/09/2011
 * @author Enrique Ramos Melgar
 * http://www.esc-studio.com
 *
 * THIS LIBRARY IS RELEASED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE
 * http://creativecommons.org/licenses/by/3.0/
 * http://www.processing.org/
 */
package BeamCalc;

import java.util.ArrayList;
import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PGraphics;
import processing.core.PVector;

/**
 * This class contains the definitions of every node in the
 * structure, the integration algorithms and the equilibrium
 * equations. Node visualisation is processed here. 
 * @author Enrique Ramos
 */
public class Node {

	
	private PApplet p; 				// Parent Processing Application
	private Structure structure; 	// Parent Structure
	
	// Elements linked to the node
	ArrayList<Element> nodeBeams = new ArrayList<Element>();
	
	// ArrayLists for the graphs
	public ArrayList<PVector> graph = new ArrayList<PVector>();
	public ArrayList<Integer> graphTimes = new ArrayList<Integer>();
	
	private float dispScale = 0;

	// Node Position Vectors
	PVector initialPos;
	public PVector pos;
	protected PVector visualPos;
	
	// Rotation vectors
	public PVector initialRot;
	public PVector rot = new PVector();
	PVector visualRot = new PVector();

	public PVector forceVec = new PVector();
	public PVector momentVec = new PVector();
	public PVector loadVec = new PVector();
	private PVector result = new PVector();
	private PVector momentResult = new PVector();
	private PVector deltaPos = new PVector();
	private PVector deltaRot = new PVector();

	PVector dispVec = new PVector();
	float rotVec;
	float windLoad;

	PVector prevPos = new PVector();
	PVector prevRot = new PVector();
	boolean switcher = false;

	public int id;
	public boolean linked = false;
	public boolean clamped = false;
	public boolean selected = false;
	private int activeGizmo = 3;

	private PVector posAcc = new PVector();
	private PVector rotAcc;
	private float posMass = 1;
	private float rotMass = 1;
	public boolean hover;

	float posStep = 0.00001f;
	private float rotStep = 0.0000000001f;
	private int startMillis;
	PVector[] damping = new PVector[2];
	PVector[] vel = new PVector[2];
	private int graphCount;

	Node(PApplet p, Structure structure, float x, float y, float z) {
		this.p = p;
		this.structure = structure;
		pos = new PVector(x, y, z);
		prevPos = new PVector(x, y, z);
		visualPos = new PVector();
		initPos();
		// loadVec.set(0, 0, -100);
	}

	Node(PApplet p, Structure structure, PVector iVec) {
		this.p = p;
		this.structure = structure;
		pos = new PVector(iVec.x, iVec.y, iVec.z);
		prevPos = new PVector(iVec.x, iVec.y, iVec.z);
		visualPos = new PVector();
		// loadVec.set(0, 0, -100);
		initPos();
		startMillis = p.millis();
		vel[0] = new PVector();
		vel[1] = new PVector();

	}


	public void euler(float h1, float h2) {

		PVector[] k1 = f(pos, rot);
		k1[0].mult(h1);
		k1[1].mult(h2);

		if (!linked)
			pos.add(k1[0]);
		if (!clamped)
			rot.add(k1[1]);

		dispVec = PVector.sub(pos, initialPos);
		graph.add(dispVec);
		graphTimes.add(p.millis() - startMillis);
	}

	public void verlet(float h1, float h2) {

		PVector[] k1 = f(pos, rot, vel, h1, h2);

		if (!linked) {
			pos.add(PVector.mult(k1[0], h1));
			// vel[0] = k1[0].get();
		}
		if (!clamped) {
			rot.add(PVector.mult(k1[1], h2));
			// vel[1] = k1[1].get();
		}

		dispVec = PVector.sub(pos, initialPos);
		if (graphCount++ % 100 == 0) {
			graph.add(dispVec);
			graphTimes.add(p.millis() - startMillis);
		}
	}


	public PVector[] f(PVector currentPos, PVector currentRot) {

		PVector posTemp = this.pos.get();
		PVector rotTemp = this.rot.get();

		PVector[] result = new PVector[2];
		result[0] = new PVector();
		result[1] = new PVector();

		forceVec = new PVector();
		momentVec = new PVector();

		pos = currentPos.get();
		rot = currentRot.get();

		// TODO Synchronous??
		for (int i = 0; i < nodeBeams.size(); i++) {
			//nodeBeams.get(i).eval();
			forceVec.add(nodeBeams.get(i).getReaction(this));
			momentVec.add(nodeBeams.get(i).getMoment(this));
		}

		// TODO
		posMass = nodeBeams.size();
		result[0] = PVector.mult(forceVec.get(), 1 / posMass);

		rotMass = nodeBeams.size();
		result[1] = PVector.mult(momentVec.get(), 1 / rotMass);

		pos = posTemp.get();
		rot = rotTemp.get();

		return result;
	}

	public PVector[] f(PVector currentPos, PVector currentRot, PVector[] cVel,
			float h1, float h2) {

		PVector posTemp = this.pos.get();
		PVector rotTemp = this.rot.get();

		PVector[] acc = new PVector[2];
		acc[0] = new PVector();
		acc[1] = new PVector();

		forceVec = new PVector();
		momentVec = new PVector();

		pos = currentPos.get();
		rot = currentRot.get();

		float mass = 0;
		float rotMass = 0;

		for (int i = 0; i < nodeBeams.size(); i++) {
			//nodeBeams.get(i).eval();
			forceVec.add(nodeBeams.get(i).getReaction(this));
			momentVec.add(nodeBeams.get(i).getMoment(this));
			mass += nodeBeams.get(i).getMass() / 2;
			rotMass += nodeBeams.get(i).getMass()
					* Math.pow(nodeBeams.get(i).L, 2) / 105;
		}

		damping[0] = PVector.mult(cVel[0], -100);
		forceVec.add(damping[0]);
		damping[1] = PVector.mult(cVel[1], -100000);
		momentVec.add(damping[1]);

		acc[0] = PVector.mult(forceVec.get(), 1 / mass);
		acc[1] = PVector.mult(momentVec.get(), 1 / mass * 100f);

		if (!linked)
			cVel[0].add(PVector.mult(acc[0], h1));
		if (!clamped)
			cVel[1].add(PVector.mult(acc[1], h2));

		pos = posTemp.get();
		rot = rotTemp.get();

		return cVel;
	}

/**
 * This is superseded by euler and verlet functions
 * @deprecated
 */
	private void getForces() {

		for (int i = 0; i < nodeBeams.size(); i++) {
			nodeBeams.get(i).eval();
		}

		forceVec = new PVector();
		momentVec = new PVector();

		for (int i = 0; i < nodeBeams.size(); i++) {
			momentVec.add(nodeBeams.get(i).getMoment(this));
			forceVec.add(nodeBeams.get(i).getReaction(this));
		}
		// Add concentrated loads and moments on nodes
		// forceVec.add(loadVec);
		// momentVec.add(concMomentVec);

	}

	/**
	 * This is superseded by euler and verlet functions
	 * @deprecated
	 */
	private void move() {

		if (linked) {
			result = PVector.mult(forceVec.get(), -1);
			forceVec.add(result);
		} else {

			posMass = nodeBeams.size();
			posAcc = PVector.mult(forceVec.get(), 1 / posMass);
			deltaPos = PVector.mult(posAcc, posStep);
			pos.add(deltaPos);

		}

	}

	/**
	 * This is superseded by euler and verlet functions
	 * @deprecated
	 */
	private void rotate() {

		if (clamped) {
			momentResult = PVector.mult(momentVec.get(), -1);
			momentVec.add(momentResult);
		} else {

			rotMass = nodeBeams.size();
			rotAcc = PVector.mult(momentVec.get(), 1 / posMass);
			deltaRot = PVector.mult(rotAcc.get(), rotStep);
			rot.add(deltaRot);
		}

	}

	/**
	 * When an element is created, it calls this function to be added
	 * in the nodeBeams arrayList in each of the nodes.
	 * @param BeamIn
	 */
	public void addBeam(Element BeamIn) {
		nodeBeams.add(BeamIn);
	}

	/**
	 * Isolates the current node
	 */
	public void unLink() {
		nodeBeams = new ArrayList<Element>();
	}
	
/**
 * Draws the current node to the PApplet passed on as a parameter
 * @param window
 */
	public void draw(PApplet window) {

		// Get the visualisation scale from the toolbox
		dispScale = structure.toolbox.getDispScale();

		// This computes the current displacement of the node from its original position
		visualPos = new PVector(initialPos.x + dispVec.x * dispScale,
				initialPos.y + dispVec.y * dispScale, initialPos.z + dispVec.z
						* dispScale);

		// This computes the current rotation of the node from its original rotation
		visualRot = new PVector(dispScale * rot.x, dispScale * rot.y, dispScale
				* rot.z);

		// Draw Rotation Lines
		window.pushMatrix();
		window.translate(visualPos.x, visualPos.y, visualPos.z);

		if (selected)
			displayData();
		window.rotateX(-visualRot.x);
		window.rotateY(-visualRot.y);
		window.rotateZ(-visualRot.z);

		// Select the colour of the node
		if (clamped) {
			window.fill(255, 0, 0);
			window.stroke(255, 0, 0);
		} else {
			if (linked) {
				window.fill(0, 0, 255);
				window.stroke(0, 0, 255);
			} else {
				window.fill(0, 255, 0);
				window.stroke(0, 255, 0);
			}
		}
		if (hover)
			window.stroke(255, 200, 0);
		
		// Draw a point
		window.strokeWeight(10);
		if(!structure.rendering)window.point(0, 0, 0);
		window.strokeWeight(1);
		window.popMatrix();

		// drawVec(forceVec, p, 1, 0);
		// drawVec(damping[0], p, 1000, 1);

		if (selected)
			drawGizmo(window.g, 500, 3, activeGizmo);
	}

	private void drawForces(PApplet window, float scale) {

		window.strokeWeight(5);
		window.stroke(255, 0, 0);
		window.line(visualPos.x, visualPos.y, visualPos.z, visualPos.x
				+ forceVec.x * scale, visualPos.y + forceVec.y * scale,
				visualPos.z + forceVec.z * scale);
		window.strokeWeight(1);

	}

	private void drawReaction(PApplet window, float scale) {

		// Draw reaction force
		window.stroke(255, 0, 255);
		window.line(visualPos.x, visualPos.y, visualPos.z, visualPos.x
				+ result.x * scale, visualPos.y + result.y * scale, visualPos.z
				+ result.z * scale);
	}

	private void drawVec(PVector vec, PApplet window, float scale, int col) {

		// Draw reaction force
		window.stroke(col * 255, 0, 0);
		window.line(visualPos.x, visualPos.y, visualPos.z, visualPos.x + vec.x
				* scale, visualPos.y + vec.y * scale, visualPos.z + vec.z
				* scale);
	}

	/**
	 * Displays the id of the node
	 */
	public void displayData() {

		p.scale(1, -1, 1);
		p.fill(0);
		p.textAlign(PConstants.LEFT);
		p.text("Node " + id, 10, 10);
		p.scale(1, -1, 1);

	}

	/**
	 * Returns the connection state of the node
	 * @return
	 */
	public boolean isDisconnected() {
		if (nodeBeams.size() < 1)
			return true;
		return false;
	}

	/**
	 * Returns the displacement vector of the node
	 * @return
	 */
	public float getDisp() {
		float disp;
		disp = dispVec.mag();
		return disp;
	}

	/**
	 * Returns the node to its initial position
	 */
	private void initPos() {
		initialPos = new PVector(pos.x, pos.y, pos.z);
	}

	/**
	 * Draws a gizmo to allow the user to modify the position and
	 * rotation of the selected node
	 * @param window
	 * @param scale
	 * @param weight
	 * @param active
	 */
	public void drawGizmo(PGraphics window, float scale, int weight, int active) {

		PVector gizmoX = new PVector(1, 0, 0);
		PVector gizmoY = new PVector(0, 1, 0);
		PVector gizmoZ = new PVector(0, 0, 1);

		if (active == 0) {
			window.stroke(230, 200, 30);
			window.strokeWeight(2 * weight);
		} else {
			window.stroke(255, 0, 0);
			window.strokeWeight(weight);
		}
		window.line(visualPos.x, visualPos.y, visualPos.z, visualPos.x
				+ gizmoX.x * scale, visualPos.y + gizmoX.y * scale, visualPos.z
				+ gizmoX.z * scale);
		window.strokeWeight(4 * weight);
		window.point(visualPos.x + gizmoX.x * scale, visualPos.y + gizmoX.y
				* scale, visualPos.z + gizmoX.z * scale);

		if (active == 1) {
			window.stroke(230, 200, 30);
			window.strokeWeight(2 * weight);
		} else {
			window.stroke(0, 255, 0);
			window.strokeWeight(weight);
		}
		window.line(visualPos.x, visualPos.y, visualPos.z, visualPos.x
				+ gizmoY.x * scale, visualPos.y + gizmoY.y * scale, visualPos.z
				+ gizmoY.z * scale);
		window.strokeWeight(4 * weight);
		window.point(visualPos.x + gizmoY.x * scale, visualPos.y + gizmoY.y
				* scale, visualPos.z + gizmoY.z * scale);

		if (active == 2) {
			window.stroke(230, 200, 30);
			window.strokeWeight(2 * weight);
		} else {
			window.stroke(0, 0, 255);
			window.strokeWeight(weight);
		}
		window.line(visualPos.x, visualPos.y, visualPos.z, visualPos.x
				+ gizmoZ.x * scale, visualPos.y + gizmoZ.y * scale, visualPos.z
				+ gizmoZ.z * scale);
		window.strokeWeight(4 * weight);
		window.point(visualPos.x + gizmoZ.x * scale, visualPos.y + gizmoZ.y
				* scale, visualPos.z + gizmoZ.z * scale);
		window.strokeWeight(1);

	}

	public void activeGizmo(int i) {
		this.activeGizmo = i;
	}

	public void modify(float dist, int mode) {

		if (mode == 0) {
			if (activeGizmo == 0)
				this.pos.x += dist / 10;
			if (activeGizmo == 1)
				this.pos.y += dist / 10;
			if (activeGizmo == 2)
				this.pos.z += dist / 10;
		} else {
			if (activeGizmo == 0)
				this.rot.x += dist / 100000;
			if (activeGizmo == 1)
				this.rot.y += dist / 100000;
			if (activeGizmo == 2)
				this.rot.z += dist / 100000;
		}

	}

	public float getYscaleGraph() {

		float maxDisp = 0;
		for (int i = 0; i < graph.size(); i++) {
			if (Math.abs(graph.get(i).z) > maxDisp)
				maxDisp = Math.abs(graph.get(i).z);
		}
		return maxDisp;
	}

	public float getXscaleGraph() {

		float maxDisp = 0;
		for (int i = 0; i < graphTimes.size(); i++) {
			if (Math.abs(graphTimes.get(i)) > maxDisp)
				maxDisp = Math.abs(graphTimes.get(i));
		}
		return maxDisp;
	}
}
