/**
 * (./) Element.java v0.1 05/09/2011
 * @author Enrique Ramos Melgar
 * http://www.esc-studio.com
 *
 * THIS LIBRARY IS RELEASED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE
 * http://creativecommons.org/licenses/by/3.0/
 * http://www.processing.org/
 */

package BeamCalc;

import processing.core.PApplet;
import processing.core.PConstants;
import processing.core.PVector;

/**
 * The class Element is where the calculation of the reactions
 * and the internal stresses of the element takes place.
 * @author Enrique Ramos
 * 
 **/
public class Element {

	private static PApplet p;			// Parent Processing Application
	private Structure structure;		// Parent Structure
	private CoordinateSystem initCS;	// Initial Local Coordinate System
	private CoordinateSystem localCS; 	// Local Coordinate System
	private CoordinateSystem visualCS;	// Coordinate System for Visualisation
	private Material mat;				// Element Material	
	Section section;					// Element Section
	Gradient elmGrad = new Gradient(); 	// Gradient
	
	// Nodes and element id
	public int id;
	public Node a;
	public Node b;
	
	public boolean selected = false; 				// Stelection State
	public PVector aAngles = new PVector(); 		// Angles of node a
	public PVector bAngles = new PVector();			// Angles of node b
	public PVector planeVec = new PVector();		// Principal Section Plane
	public PVector beamAngles = new PVector();		// Angles of the element
	public PVector beamInitAngles = new PVector();	// Initial Angles of the Element
	
	PVector beamVec = new PVector();				// Vector defining the element
	PVector beamVecLocal = new PVector();			// Vector defining the element in local CS
	PVector weightVec = new PVector();				// Vector defining the weight of the element

	public PVector aMoment = new PVector();			// Member-end moments at node a
	public PVector bMoment = new PVector();			// Member-end moments at node b
	public PVector aReact = new PVector();			// Member-end forces at node a
	public PVector bReact = new PVector();			// Member-end forces at node b

	// Section Properties
	private float E; // modulus of elasticity (N/mm2)
	private float Ixz; // moment of inertia (mm4)
	private float Ixy; // moment of inertia (mm4)
	private float J; // moment of inertia in Torsion (mm4)
	private float G; // Shear Modulus (N/mm2)
	private float A; // Area of Section (mm2)
	private float D; // Density (T/mm2)

	public float L; // Span length (mm)
	public float initialL; // Length of bar when relaxed (mm)
	private float delta; // Length increase (mm)
	private float deltaRot; // Difference in Rotation (Rad)

	public float lW; // Element Linear Weight (N/m)
	public float W; // Element total Weight

	int numPoints; // number of points at which to compute values

	int i; // Variable for loops

	// arrays to store distance, shear, moment and deflection
	float x[]; // distance from the left end of the beam
	float v[][]; // shear force
	float m[][]; // moment
	float y[][]; // deflection
	float d[][]; // slope
	float rr1[]; // reaction at left support
	float rr2[]; // reaction at right support
	float mm1[]; // moment at left end
	float mm2[]; // moment at right end

	public float axial;		// Axial Stress
	float weightAxial; 		// Axial component of weight
	public float torsion;	// Element Torsion

	// arrays to store distance, shear, moment and deflection from Weight
	float wv[][]; // shear force
	float wm[][]; // moment
	float wy[][]; // deflection
	float wd[][]; // slope
	float wrr1[]; // reaction at left support
	float wrr2[]; // reaction at right support
	float wmm1[]; // moment at left end
	float wmm2[]; // moment at right end

	public boolean fullCalc;
	float dispLScale;
	public boolean render = false;
	private int[] elmCol = new int[3];
	private float elmStroke;
	public boolean hover;
	public boolean clamped;

	/**
	 *  Main constructor of the class Element
	 * @param pa: Parent Processing PApplet
	 * @param structure: Parent Structure Instance
	 * @param a: Node a
	 * @param b: Node b
	 * @param numPoints: Number of point for internal stress analysis
	 * @param planeVec: Main Section Plane
	 */
	public Element(PApplet pa, Structure structure, Node a, Node b,
			int numPoints, PVector planeVec) {

		Element.p = pa;
		this.structure = structure;
		this.a = a;
		this.b = b;

		this.planeVec = planeVec;

		// Add this beam to the node beam array
		a.addBeam(this);
		b.addBeam(this);

		// Initialise Coordinate Systems
		initCS = new CoordinateSystem();
		localCS = new CoordinateSystem();
		visualCS = new CoordinateSystem();

		beamVec = PVector.sub(b.pos, a.pos);
		initCS.alignXToVector(beamVec, planeVec);
		localCS.alignXToVector(beamVec, planeVec);

		beamVecLocal = initCS.changeTo(beamVec);

		localCS.pos = a.pos.get();

		aAngles = localCS.changeTo(a.rot);
		bAngles = localCS.changeTo(b.rot);

		this.initialL = PVector.dist(a.pos, b.pos);
		this.L = initialL;

		// Material Properties
		this.mat = new Material();
		this.E = mat.getE(); // N/mm2
		this.G = mat.getG(); // N/mm2
		this.D = mat.getD(); // Density (N/mm3)

		// Section Properties
		this.section = new Section();
		this.Ixz = section.getIxz(); // mm4 1338000
		this.Ixy = section.getIxy(); // mm4 1338000
		this.J = section.getJ(); // mm4
		this.A = section.getA(); // mm2

		this.lW = A * D;
		this.numPoints = numPoints;

		x = new float[numPoints];

		v = new float[2][numPoints];
		m = new float[2][numPoints];
		y = new float[2][numPoints];
		d = new float[2][numPoints];

		wv = new float[2][numPoints];
		wm = new float[2][numPoints];
		wy = new float[2][numPoints];
		wd = new float[2][numPoints];

		rr1 = new float[2];
		rr2 = new float[2];
		mm1 = new float[2];
		mm2 = new float[2];

		wrr1 = new float[2];
		wrr2 = new float[2];
		wmm1 = new float[2];
		wmm2 = new float[2];

		// Set the points that will be calculated
		for (i = 0; i < numPoints; ++i)
			x[i] = i * L / ((float) numPoints - 1);

		evalWeight();
	}

	/**
	 * Evaluation of the internal Stresses and reactions
	 * of the element
	 */
	public void eval() {

		this.Ixz = section.getIxz(); // mm4 
		this.Ixy = section.getIxy(); // mm4 
		
		// Calculate the angle of the bar
		beamVec = PVector.sub(b.pos, a.pos);
		beamVecLocal = initCS.changeTo(beamVec);
		L = PVector.dist(a.pos, b.pos);

		beamAngles.z = -(float) (Math.atan2(beamVecLocal.y, beamVecLocal.x));
		beamAngles.y = -(float) (Math.atan2(beamVecLocal.z, beamVecLocal.x));

		localCS.alignXToVector(beamVec, planeVec);
		localCS.pos = a.pos.get();

		// Calculate Node Angles in Local Coordinates
		aAngles = localCS.changeTo(a.rot);
		bAngles = localCS.changeTo(b.rot);
		bAngles.y *= -1;
		aAngles.y *= -1;
		bAngles.x *= -1;
		aAngles.x *= -1;

		// Scale the display length according to the nes Length
		dispLScale = PVector.dist(a.visualPos, b.visualPos) / L;

		// Maybe this could be left out for small deformations?
		for (i = 0; i < numPoints; ++i)
			x[i] = i * L / ((float) numPoints - 1);

		// These are the main calculations
		evalAxial();	
		evalRot();
		evalMomentXZ();
		evalMomentXY();

		// This calculates the final member-end forces
		aMoment.x = torsion;
		aMoment.y = -mm1[0];
		aMoment.z = mm1[1];

		bMoment.x = -torsion;
		bMoment.y = mm2[0];
		bMoment.z = -mm2[1];

		aReact.x = axial - weightAxial;
		aReact.y = rr1[1];
		aReact.z = rr1[0];

		bReact.x = -axial - weightAxial;
		bReact.y = rr2[1];
		bReact.z = rr2[0];

	}

	private void evalMomentXZ() {
		evalFromAngles(aAngles.y - beamAngles.y, -bAngles.y + beamAngles.y,
				Ixz, 0);
	}

	private void evalMomentXY() {
		evalFromAngles(aAngles.z - beamAngles.z, -bAngles.z + beamAngles.z,
				Ixy, 1);
	}

	private void evalFromAngles(float slopeA, float slopeB, float I, int axis) {

		float ma = -(2 * slopeA - slopeB) * 2 * E * I / L;
		float mb = (2 * slopeB - slopeA) * 2 * E * I / L;

		float rr1Temp, rr2Temp, mm1Temp, mm2Temp;

		reactMoment(ma, 0, axis);

		rr1Temp = rr1[axis];
		rr2Temp = rr2[axis];
		mm1Temp = mm1[axis];
		mm2Temp = mm2[axis];

		reactMoment(mb, L, axis);
		rr1[axis] += rr1Temp;
		rr2[axis] += rr2Temp;
		mm1[axis] += mm1Temp;
		mm2[axis] += mm2Temp;

		addWeightReact(axis);

	}

	private void reactMoment(float mo, float a, int axis) {
		// compute reaction
		rr1[axis] = -mo / L;
		rr2[axis] = -rr1[axis];
		if (a == L) {
			mm1[axis] = 0;
			mm2[axis] = rr1[axis] * L;
		} else {
			mm1[axis] = mo;
			mm2[axis] = rr1[axis] * L + mo;
		}

	}

	/**
	 * This Method evaluates internal stresses caused by a moment applied at a
	 * distance a from the starting node, on the selected axis.
	 */
	private void evalMoment(float mo, float a, float I, int axis) {

		for (i = 0; i < numPoints; ++i) {

			if (a == L) {
				v[axis][i] = rr1[axis];
				m[axis][i] = rr1[axis] * x[i];
				y[axis][i] = (float) (mo / (6 * E * L * I) * ((6 * a * L - 3
						* a * a - 2 * L * L)
						* x[i] - Math.pow(x[i], 3)));
				d[axis][i] = (float) (mo / (6 * E * L * I) * ((6 * a * L - 3
						* a * a - 2 * L * L) - 3 * Math.pow(x[i], 2)));
			} else {
				v[axis][i] = rr1[axis];
				m[axis][i] = rr1[axis] * x[i] + mo;
				y[axis][i] = (float) (mo / (6 * E * L * I) * (3 * a * a * L + 3
						* Math.pow(x[i], 2) * L - Math.pow(x[i], 3) - (2 * L
						* L + 3 * a * a)
						* x[i]));
				d[axis][i] = (float) (mo / (6 * E * L * I) * (6 * x[i] * L - 3
						* Math.pow(x[i], 2) - (2 * L * L + 3 * a * a)));
			}

		}
	}

	void evalWeight() {

		weightVec = localCS.changeTo(PVector.mult(
				CoordinateSystem.universalCS.zAxis, lW));

		weightAxial = weightVec.x * L / 2;

		reactLinear(weightVec.z, 0);
		evalLinear(weightVec.z, Ixz, 0);

		reactLinear(weightVec.y, 1);
		evalLinear(weightVec.y, Ixy, 1);

		// float vertAngle = PVector.angleBetween(beamVec,
		// CoordinateSystem.universalCS.zAxis);
		// float weightNormal = (float) Math.abs((lW * Math.sin(vertAngle)));
		// weightAxial = (float) -(lW * Math.cos(vertAngle) * L / 2);
		// reactLinear(weightNormal, 0);
		// evalLinear(weightNormal, Ixz, 0);
	}

	private void addWeightReact(int axis) {

		rr1[axis] += wrr1[axis];
		rr2[axis] += wrr2[axis];
		mm1[axis] += wmm1[axis];
		mm2[axis] += wmm2[axis];
	}

	private void addWeight(int axis) {
		for (i = 0; i < numPoints; ++i) {
			v[axis][i] += wv[axis][i];
			m[axis][i] += wm[axis][i];
			y[axis][i] += wy[axis][i];
			d[axis][i] += wd[axis][i];
		}
	}

	public void updateInternal() {

		updateElmData(aAngles.y - beamAngles.y, -bAngles.y + beamAngles.y, Ixz,
				0);
		updateElmData(aAngles.z - beamAngles.z, -bAngles.z + beamAngles.z, Ixy,
				1);

	}

	private void updateElmData(float slopeA, float slopeB, float I, int axis) {

		float ma = -(2 * slopeA - slopeB) * 2 * E * I / L;
		float mb = (2 * slopeB - slopeA) * 2 * E * I / L;

		float[] vTemp = new float[numPoints];
		float[] mTemp = new float[numPoints];
		float[] yTemp = new float[numPoints];
		float[] dTemp = new float[numPoints];
		float rr1Temp, rr2Temp, mm1Temp, mm2Temp;

		reactMoment(ma, 0, axis);
		evalMoment(ma, 0, I, axis);
		rr1Temp = rr1[axis];
		rr2Temp = rr2[axis];
		mm1Temp = mm1[axis];
		mm2Temp = mm2[axis];
		for (i = 0; i < numPoints; ++i) {
			vTemp[i] = v[axis][i];
			mTemp[i] = m[axis][i];
			yTemp[i] = y[axis][i];
			dTemp[i] = d[axis][i];
		}
		reactMoment(mb, L, axis);
		evalMoment(mb, L, I, axis);
		rr1[axis] += rr1Temp;
		rr2[axis] += rr2Temp;
		mm1[axis] += mm1Temp;
		mm2[axis] += mm2Temp;
		for (i = 0; i < numPoints; ++i) {
			v[axis][i] += vTemp[i];
			m[axis][i] += mTemp[i];
			y[axis][i] += yTemp[i];
			d[axis][i] += dTemp[i];
		}

		addWeight(axis);

	}

	private void reactLinear(float w1, int axis) {

		W = -w1 * L; // total load
		wrr1[axis] = W / 2; // reaction
		wrr2[axis] = wrr1[axis];
		wmm1[axis] = -W * L / 12;
		wmm2[axis] = wmm1[axis]; // reaction
	}

	private void evalLinear(float w1, float I, int axis) {

		for (i = 0; i < numPoints; ++i) {

			wv[axis][i] = wrr1[axis] - W * x[i] / L;

			wm[axis][i] = (float) (wmm1[axis] + wrr1[axis] * x[i] - 0.5 * W
					* x[i] * x[i] / L);

			wy[axis][i] = (float) (1. / (6. * E * I) * (wrr1[axis]
					* Math.pow(x[i], 3.) - 3. * -wmm1[axis] * x[i] * x[i] - W
					/ 4. * Math.pow(x[i], 4) / L));

			wd[axis][i] = (float) (1. / (6. * E * I) * (wrr1[axis] * 3
					* Math.pow(x[i], 2.) - 6. * -wmm1[axis] * x[i] - W / 4. * 4
					* Math.pow(x[i], 3.) / L));
		}
	}

	/**
	 * This Method evaluates the axial stress Tension is positive, Compression
	 * Negative
	 */
	private void evalAxial() {

		delta = L - initialL;

		axial = E * A * delta / initialL;
	}

	/**
	 * This Method evaluates the rotational stress
	 */
	private void evalRot() {

		deltaRot = aAngles.x - bAngles.x;
		torsion = deltaRot * G * J / L;
	}

	/**
	 * Get the member-end moment vector at node "node"
	 * @param node
	 * @return
	 */
	public PVector getMoment(Node node) {
		// p.println("a"+this.mm1 + "   b"+this.mm2);

		if (node == a) {
			return localCS.changeFrom(aMoment);
		} else {
			return localCS.changeFrom(bMoment);
		}
	}

	/**
	 * Get the member-end force at node "node"
	 */
	public PVector getReaction(Node node) {
		if (node == a)
			return this.localCS.changeFrom(aReact);
		if (node == b)
			return this.localCS.changeFrom(bReact);
		return null;
	}
	
	/**
	 * Get the maximum moment value in the element
	 * @return
	 */
	public float getMaxMoment() {

		float maxMoment = 0;
		for (int j = 0; j < 2; j++) {
			for (int i = 0; i < numPoints; i++) {
				if (Math.abs(m[j][i]) > maxMoment)
					maxMoment = Math.abs(m[j][i]);
			}
		}
		return maxMoment;
	}
	
	/**
	 * Get maximum shear value in the element
	 * @return
	 */
	public float getMaxShear() {

		float maxShear = 0;
		for (int j = 0; j < 2; j++) {
			for (int i = 0; i < numPoints; i++) {
				if (Math.abs(v[j][i]) > maxShear)
					maxShear = Math.abs(v[j][i]);
			}
		}
		return maxShear;
	}

	/**
	 * Get maximum node displacement value
	 * @return
	 */
	public float getMaxNodeDisp() {
		if (a.dispVec.mag() > b.dispVec.mag()) {
			return a.dispVec.mag();
		} else {
			return b.dispVec.mag();
		}
	}
	
	/**
	 * Get the absolute maximum element displacement value
	 * @return
	 */
	public float getMaxDisp() {

		float maxDisp = 0;
		for (int j = 0; j < 2; j++) {
			for (int i = 0; i < numPoints; i++) {
				if (Math.abs(y[j][i]) > maxDisp)
					maxDisp = Math.abs(y[j][i]);
			}
		}
		return maxDisp;
	}
	
	/**
	 * Returns the element's mass
	 * @return
	 */
	public double getMass(){
		return lW*initialL*0.1f;
	}

	/**
	 * Returns the opposite node to the one passed on
	 * @param thisNode
	 * @return
	 */
	public Node getOtherNode(Node thisNode) {
		if (thisNode == a)
			return b;
		return a;
	}

	/**
	 * Draw the element
	 */
	public void draw() {

		visualCS.alignXToVector(PVector.sub(b.visualPos, a.visualPos), planeVec);
		visualCS.pos = a.visualPos.get();

		if (selected) {
			p.strokeWeight(3);
			visualCS.drawCS(p, 1000);
			p.strokeWeight(1);
		}

		// Draw from the beam Coordinate System
		p.pushMatrix();
		visualCS.alignMatrix(p);
		p.rotateX(aAngles.x * structure.toolbox.getDispScale());
		drawDisp(-structure.toolbox.getDispScale());

		if (selected) {
			displayData();
		}

		//drawReactions(0.01f);
		p.popMatrix();

	}

	/**
	 * Display element data
	 */
	public void displayData() {

		drawMoment(structure.toolbox.getMomentScale() / 1000);
		drawTorsion(structure.toolbox.getMomentScale() / 1000);

		p.fill(100);
		p.textAlign(PConstants.CENTER);
		p.text("Element " + id, L / 2, 20);

	}

	private void drawShear(float scale) {
		p.stroke(0, 255, 0);
		p.fill(100, 50);
		// p.noFill();
		p.beginShape();
		p.vertex(0, 0);
		for (int i = 0; i < numPoints; i++) {
			p.vertex(x[i] * dispLScale, v[1][i] * -scale, 0);
		}
		p.vertex(L * dispLScale, 0);
		p.endShape(PConstants.CLOSE);

		p.stroke(255, 0, 0);
		p.fill(100, 50);
		// p.noFill();
		p.beginShape();
		p.vertex(0, 0);
		for (int i = 0; i < numPoints; i++) {
			p.vertex(x[i] * dispLScale, 0, v[0][i] * -scale);
		}
		p.vertex(L * dispLScale, 0);
		p.endShape(PConstants.CLOSE);

	}

	private void drawMomentVec(float scale) {

		PVector zero = new PVector();
		PVector bVec = new PVector(dispLScale * L, 0, 0);
		p.stroke(200, 2, 150);
		drawVec(aMoment, zero, scale, 4);
		p.stroke(0, 250, 80);
		drawVec(bMoment, bVec, scale, 4);

	}

	private void drawReactions(float scale) {

		PVector zero = new PVector();
		PVector bVec = new PVector(dispLScale * L, 0, 0);
		p.stroke(200, 2, 150);
		drawVec(aReact, zero, scale, 4);
		p.stroke(0, 250, 80);
		drawVec(bReact, bVec, scale, 4);

	}

	public void drawVec(PVector myVec, PVector myPos, float scale, float weight) {

		p.strokeWeight(weight);
		p.line(myPos.x, myPos.y, myPos.z, myPos.x + myVec.x * scale, myPos.y
				+ myVec.y * scale, myPos.z + myVec.z * scale);
		p.strokeWeight(1);
	}

	private void drawDisp(float scale) {

		p.noFill();
		if (selected) {
			p.stroke(0, 0, 0);
			p.strokeWeight(4);
		} else {
			p.stroke(0);
			p.strokeWeight(2);
		}

		p.fill(255);

		if (structure.displayMode != 0) {
			p.strokeWeight(elmStroke);
			p.fill(elmCol[0], elmCol[1], elmCol[2]);
			p.stroke(elmCol[0], elmCol[1], elmCol[2]);
		}

		if (hover) {
			p.stroke(220, 200, 30);
			p.strokeWeight(4);
		}

		if (render) {
			if(hover)p.fill(200,150,10);
			if (structure.displayMode != 0){
				p.noStroke();
			}else{
				//p.strokeWeight(0.5f);
				p.stroke(50);
			}
			section.render(p, this, scale);

		} else {
			for (int i = 0; i < numPoints - 1; i++) {
				p.line(x[i] * dispLScale, y[1][i] * scale, y[0][i] * scale,
						x[i + 1] * dispLScale, y[1][i + 1] * scale, y[0][i + 1]
								* scale);
			}
		}

		p.strokeWeight(1);

	}

	private void drawMoment(float scale) {

		p.stroke(0, 255, 0);
		p.fill(100, 50);
		// p.noFill();
		p.beginShape(PConstants.QUAD_STRIP);
		for (int i = 0; i < numPoints; i++) {
			p.vertex(x[i] * dispLScale, 0, 0);
			p.vertex(x[i] * dispLScale, m[1][i] * scale, 0);
		}
		p.endShape(PConstants.CLOSE);

		p.stroke(255, 0, 0);
		// p.noStroke();
		p.fill(100, 50);
		// p.noFill();
		p.beginShape(PConstants.QUAD_STRIP);
		for (int i = 0; i < numPoints; i++) {
			p.vertex(x[i] * dispLScale, 0, 0);
			p.vertex(x[i] * dispLScale, 0, m[0][i] * scale);
		}
		p.endShape(PConstants.CLOSE);
	}

	private void drawTorsion(float scale) {

		p.stroke(0, 0, 255);
		p.fill(100, 50);
		p.rect(0, 0, L * dispLScale, torsion * scale);

	}

	private void drawSlopes(float scale) {
		//
		// // Show Tangent
		// float lineSize = 100;
		// p.stroke(0, 255, 0);
		// for (int i = 0; i < numPoints; i += 2) {
		// p.line(x[i] * dispLScale, y[i] * scale, x[i] * dispLScale
		// + lineSize * (float) Math.cos(d[i]), y[i] * scale
		// + lineSize * (float) Math.sin(d[i]) * scale);
		// }
	}

	public PVector trim(Node node, PVector deltaPos) {

		PVector trimmed = localCS.changeTo(deltaPos);
		if (node == a) {
			float cDist = b.pos.dist(PVector.add(a.pos, deltaPos));
			trimmed.x -= (initialL - cDist) * 0.99;
		} else {
			float cDist = a.pos.dist(PVector.add(b.pos, deltaPos));
			trimmed.x += (initialL - cDist) * 0.99;
		}
		return localCS.changeFrom(trimmed);
	}

	public void setColorRange(int mode, float min, float max) {

		// TODO Gradient properly in three colours

		float t = 0;

		switch (mode) {
		case 0:
			break;
		case 1:
			t = axial;
			if (Math.abs(structure.maxAxial) > Math.abs(structure.minAxial)) {
				elmStroke = 2 + 5 * Math.abs(t / structure.maxAxial);
			} else {
				elmStroke = 2 + 5 * Math.abs(t / structure.minAxial);
			}
			break;
		case 2:
			t = getMaxMoment();
			elmStroke = 2 + 5 * (t - structure.minMoment)
					/ (structure.maxMoment - structure.minMoment);
			break;
		case 3:
			t = Math.abs(torsion);
			elmStroke = 2 + 5 * Math.abs((t - structure.minTorsion)
					/ (structure.maxTorsion - structure.minTorsion));
			break;
		case 4:
			t = getMaxNodeDisp();
			elmStroke = 2 + 5 * Math.abs((t - structure.minDisp)
					/ (structure.maxDisp - structure.minDisp));
			break;

		}

		elmGrad.setBounds(min, max);

		elmCol = elmGrad.getColor(t);

	}

	public void rotateX(int k) {
		initCS.rotateX(k);
		
	}

	public void reset() {
		
		x = new float[numPoints];

		v = new float[2][numPoints];
		m = new float[2][numPoints];
		y = new float[2][numPoints];
		d = new float[2][numPoints];

		wv = new float[2][numPoints];
		wm = new float[2][numPoints];
		wy = new float[2][numPoints];
		wd = new float[2][numPoints];

		rr1 = new float[2];
		rr2 = new float[2];
		mm1 = new float[2];
		mm2 = new float[2];

		wrr1 = new float[2];
		wrr2 = new float[2];
		wmm1 = new float[2];
		wmm2 = new float[2];

		// Set the points that will be calculated
		for (i = 0; i < numPoints; ++i)
			x[i] = i * L / ((float) numPoints - 1);
		
	}
}
