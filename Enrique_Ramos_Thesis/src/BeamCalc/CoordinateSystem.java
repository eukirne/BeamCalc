/**
 * (./) CoordinateSystem.java v0.1 05/09/2011
 * @author Enrique Ramos Melgar
 * http://www.esc-studio.com
 *
 * THIS LIBRARY IS RELEASED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE
 * http://creativecommons.org/licenses/by/3.0/
 * http://www.processing.org/
 */

package BeamCalc;

import processing.core.PApplet;
import processing.core.PMatrix3D;
import processing.core.PVector;

/**
 * This class implements a Coordinate System, and contains
 * the methods necessary for the translation of one vector
 * from a CS to another
 * @author Enrique Ramos
 *
 */
public class CoordinateSystem {

	public static CoordinateSystem universalCS = new CoordinateSystem(
			new PVector(1, 0, 0), new PVector(0, 1, 0), new PVector(0, 0, 1));

	public PVector xAxis;
	public PVector yAxis;
	public PVector zAxis;

	private float[][] CSArray = new float[4][4];
	private PMatrix3D CSMatrix;

	public PVector pos = new PVector();

	public CoordinateSystem() {

		xAxis = new PVector(1, 0, 0);
		yAxis = new PVector(0, 1, 0);
		zAxis = new PVector(0, 0, 1);

	}

	/**
	 * Constructor from three vectors
	 * @param xAxis
	 * @param yAxis
	 * @param zAxis
	 */
	public CoordinateSystem(PVector xAxis, PVector yAxis, PVector zAxis) {

		this.xAxis = xAxis;
		this.yAxis = yAxis;
		this.zAxis = zAxis;
	}

	/**
	 * Align the x axis of the CS to the input vector
	 * @param inputVec
	 */
	public void alignXToVector(PVector inputVec) {

		xAxis = inputVec.get();
		if (xAxis.x == 0 && xAxis.y == 0) {
			yAxis = universalCS.yAxis.get();
		} else {
			yAxis = universalCS.zAxis.cross(xAxis);
		}
		zAxis = xAxis.cross(yAxis);
		xAxis.normalize();
		yAxis.normalize();
		zAxis.normalize();
	}

	/**
	 * Align the CS to the first parameter vector, and the Y axis to
	 * the second parameter vector
	 * @param inputVec
	 * @param planeVec
	 */
	public void alignXToVector(PVector inputVec, PVector planeVec) {

		xAxis = inputVec.get();
		zAxis = xAxis.cross(planeVec);
		yAxis = zAxis.cross(xAxis);
		xAxis.normalize();
		yAxis.normalize();
		zAxis.normalize();
	}

	/**
	 * Translates the input vector from the universal CS, to the current CS
	 * 
	 */
	public PVector changeTo(PVector globalVec) {

		PVector localVec = changeCS(globalVec, universalCS, this);
		return localVec;
	}

	/**
	 * Translates the input vector from this CS to the universal CS
	 * @param localVec
	 * @return
	 */
	public PVector changeFrom(PVector localVec) {

		PVector globalVec = changeCS(localVec, this, universalCS);
		return globalVec;
	}

	/**
	 * Translate the input vector from the local CS to global CS
	 * @param localVec
	 * @param localCS
	 * @param globalCS
	 * @return
	 */
	public PVector changeCS(PVector localVec, CoordinateSystem localCS,
			CoordinateSystem globalCS) {

		PVector newVec = new PVector();

		newVec.x = (float) (Math.cos(PVector.angleBetween(globalCS.xAxis,
				localCS.xAxis))
				* localVec.x
				+ Math.cos(PVector.angleBetween(globalCS.xAxis, localCS.yAxis))
				* localVec.y + Math.cos(PVector.angleBetween(globalCS.xAxis,
				localCS.zAxis)) * localVec.z);
		newVec.y = (float) (Math.cos(PVector.angleBetween(globalCS.yAxis,
				localCS.xAxis))
				* localVec.x
				+ Math.cos(PVector.angleBetween(globalCS.yAxis, localCS.yAxis))
				* localVec.y + Math.cos(PVector.angleBetween(globalCS.yAxis,
				localCS.zAxis)) * localVec.z);
		newVec.z = (float) (Math.cos(PVector.angleBetween(globalCS.zAxis,
				localCS.xAxis))
				* localVec.x
				+ Math.cos(PVector.angleBetween(globalCS.zAxis, localCS.yAxis))
				* localVec.y + Math.cos(PVector.angleBetween(globalCS.zAxis,
				localCS.zAxis)) * localVec.z);

		return newVec;
	}

	/**
	 * Aligns the current matrix of the PApplet p to this CS
	 * @param p
	 */
	public void alignMatrix(PApplet p) {

		// TODO This might be too slow!
		updateMatrix();
		p.translate(pos.x, pos.y, pos.z);
		p.applyMatrix(CSMatrix);
	}

	/**
	 * Translates the CS to the position newPos
	 * @param newPos
	 */
	public void translateCS(PVector newPos) {
		this.pos = newPos.get();
	}

	public void alignToCS(PApplet p) {

		// Align the x axis to the direction of the bar, and draw it
		PVector dirVec = xAxis.get();
		PVector new_up = universalCS.xAxis.get();
		new_up.normalize();
		PVector crss = dirVec.cross(new_up);

		float theAngle = PVector.angleBetween(dirVec, new_up);
		crss.normalize();

		// p.translate(pos.x,pos.y,pos.z);
		p.rotate(-theAngle, crss.x, crss.y, crss.z);
	}

	void drawCS(PApplet p, float csScale) {

		p.pushMatrix();

		p.translate(pos.x, pos.y, pos.z);
		p.scale(csScale, csScale, csScale);

		p.stroke(255, 0, 0);
		p.line(0, 0, 0, xAxis.x, xAxis.y, xAxis.z);
		p.stroke(0, 255, 0);
		p.line(0, 0, 0, yAxis.x, yAxis.y, yAxis.z);
		p.stroke(0, 0, 255);
		p.line(0, 0, 0, zAxis.x, zAxis.y, zAxis.z);
		p.popMatrix();
	}

	/**
	 * Creates a transformation matrix based on this CS
	 */
	private void updateMatrix() {

		CoordinateSystem globalCS = universalCS;
		CoordinateSystem localCS = this;

		CSArray[0][0] = (float) Math.cos(PVector.angleBetween(globalCS.xAxis,
				localCS.xAxis));
		CSArray[0][1] = (float) Math.cos(PVector.angleBetween(globalCS.xAxis,
				localCS.yAxis));
		CSArray[0][2] = (float) Math.cos(PVector.angleBetween(globalCS.xAxis,
				localCS.zAxis));
		CSArray[0][3] = 0;

		CSArray[1][0] = (float) Math.cos(PVector.angleBetween(globalCS.yAxis,
				localCS.xAxis));
		CSArray[1][1] = (float) Math.cos(PVector.angleBetween(globalCS.yAxis,
				localCS.yAxis));
		CSArray[1][2] = (float) Math.cos(PVector.angleBetween(globalCS.yAxis,
				localCS.zAxis));
		CSArray[1][3] = 0;

		CSArray[2][0] = (float) Math.cos(PVector.angleBetween(globalCS.zAxis,
				localCS.xAxis));
		CSArray[2][1] = (float) Math.cos(PVector.angleBetween(globalCS.zAxis,
				localCS.yAxis));
		CSArray[2][2] = (float) Math.cos(PVector.angleBetween(globalCS.zAxis,
				localCS.zAxis));
		CSArray[2][3] = 0;

		CSArray[3][0] = 0;
		CSArray[3][1] = 0;
		CSArray[3][2] = 0;
		CSArray[3][3] = 1;

		CSMatrix = new PMatrix3D(CSArray[0][0], CSArray[0][1], CSArray[0][2],
				CSArray[0][3], CSArray[1][0], CSArray[1][1], CSArray[1][2],
				CSArray[1][3], CSArray[2][0], CSArray[2][1], CSArray[2][2],
				CSArray[2][3], CSArray[3][0], CSArray[3][1], CSArray[3][2],
				CSArray[03][3]);

	}

	public void rotateX(int k) {
		
		PVector newVector = new PVector(0,(float)Math.cos(PApplet.radians(k)),(float) Math.sin(PApplet.radians(k)));
		this.alignXToVector(xAxis, newVector);
		
	}
}
