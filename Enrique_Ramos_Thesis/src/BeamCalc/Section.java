/**
 * (./) Section.java v0.1 05/09/2011
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

/**
 * This class contains the geometric data of the section
 * @author Enrique Ramos
 *
 */
public class Section {

	// Section Properties
	public float Ixz; // moment of inertia on the main axis (mm4)
	public float Ixy; // moment of inertia on the secondary axis(mm4)
	public float J; // moment of inertia in Torsion (mm4)
	public float A; // Area of Section (mm2)

	float beamHeight = 150;
	float beamWidth = 75;

	public Section() {

		this.Ixz = 3492000; // mm4 1338000
		// this.Ixy = 3492000; // mm4 1338000
		this.Ixy = 1338000; // mm4 1338000
		this.J = 52370000; // mm4
		this.A = 2124; // mm2

	}

	/**
	 * This function is called from the element's draw() function, and renders the current
	 * section accurately, following the element's deformation state.
	 * @param p
	 * @param elm
	 * @param scale
	 */
	public void render(PApplet p, Element elm, float scale) {

		float tempAngleXY = 0;
		float tempAngleXZ = 0;
		float tempAngleX = 0;

		p.beginShape(PConstants.QUAD_STRIP);
		for (int i = 0; i < elm.numPoints; i++) {
			if (i < elm.numPoints - 1) {
				tempAngleXZ = (float) (Math.atan2(elm.y[0][i + 1] * scale
						- elm.y[0][i] * scale, elm.x[i + 1] * elm.dispLScale
						- elm.x[i] * elm.dispLScale));
				tempAngleXY = (float) (Math.atan2(elm.y[1][i + 1] * scale
						- elm.y[1][i] * scale, elm.x[i + 1] * elm.dispLScale
						- elm.x[i] * elm.dispLScale));
				// TODO Include Torsion!!
				tempAngleX = (elm.bAngles.x * elm.L / elm.x[i] - elm.aAngles.x);
			}

			p.vertex(
					elm.x[i] * elm.dispLScale + beamHeight
							* PApplet.sin(tempAngleXZ) + beamWidth
							* PApplet.sin(tempAngleXY), elm.y[1][i] * scale
							- beamWidth * PApplet.cos(tempAngleXY), elm.y[0][i]
							* scale - beamHeight * PApplet.cos(tempAngleXZ));

			p.vertex(
					elm.x[i] * elm.dispLScale - beamHeight
							* PApplet.sin(tempAngleXZ) + beamWidth
							* PApplet.sin(tempAngleXY), elm.y[1][i] * scale
							- beamWidth * PApplet.cos(tempAngleXY), elm.y[0][i]
							* scale + beamHeight * PApplet.cos(tempAngleXZ));

		}
		p.endShape(PConstants.CLOSE);

		p.beginShape(PConstants.QUAD_STRIP);
		for (int i = 0; i < elm.numPoints; i++) {
			if (i < elm.numPoints - 1) {
				tempAngleXZ = (float) (Math.atan2(elm.y[0][i + 1] * scale
						- elm.y[0][i] * scale, elm.x[i + 1] * elm.dispLScale
						- elm.x[i] * elm.dispLScale));
				tempAngleXY = (float) (Math.atan2(elm.y[1][i + 1] * scale
						- elm.y[1][i] * scale, elm.x[i + 1] * elm.dispLScale
						- elm.x[i] * elm.dispLScale));
			}
			p.vertex(
					elm.x[i] * elm.dispLScale + beamHeight
							* PApplet.sin(tempAngleXZ) - beamWidth
							* PApplet.sin(tempAngleXY), elm.y[1][i] * scale
							+ beamWidth * PApplet.cos(tempAngleXY), elm.y[0][i]
							* scale - beamHeight * PApplet.cos(tempAngleXZ));
			p.vertex(
					elm.x[i] * elm.dispLScale - beamHeight
							* PApplet.sin(tempAngleXZ) - beamWidth
							* PApplet.sin(tempAngleXY), elm.y[1][i] * scale
							+ beamWidth * PApplet.cos(tempAngleXY), elm.y[0][i]
							* scale + beamHeight * PApplet.cos(tempAngleXZ));

		}
		p.endShape(PConstants.CLOSE);

		p.beginShape(PConstants.QUAD_STRIP);
		for (int i = 0; i < elm.numPoints; i++) {
			if (i < elm.numPoints - 1) {
				tempAngleXZ = (float) (Math.atan2(elm.y[0][i + 1] * scale
						- elm.y[0][i] * scale, elm.x[i + 1] * elm.dispLScale
						- elm.x[i] * elm.dispLScale));
				tempAngleXY = (float) (Math.atan2(elm.y[1][i + 1] * scale
						- elm.y[1][i] * scale, elm.x[i + 1] * elm.dispLScale
						- elm.x[i] * elm.dispLScale));

			}
			p.vertex(
					elm.x[i] * elm.dispLScale + beamHeight
							* PApplet.sin(tempAngleXZ) + beamWidth
							* PApplet.sin(tempAngleXY), elm.y[1][i] * scale
							- beamWidth * PApplet.cos(tempAngleXY), elm.y[0][i]
							* scale - beamHeight * PApplet.cos(tempAngleXZ));
			p.vertex(
					elm.x[i] * elm.dispLScale + beamHeight
							* PApplet.sin(tempAngleXZ) - beamWidth
							* PApplet.sin(tempAngleXY), elm.y[1][i] * scale
							+ beamWidth * PApplet.cos(tempAngleXY), elm.y[0][i]
							* scale - beamHeight * PApplet.cos(tempAngleXZ));

		}
		p.endShape(PConstants.CLOSE);

		p.beginShape(PConstants.QUAD_STRIP);
		for (int i = 0; i < elm.numPoints; i++) {
			if (i < elm.numPoints - 1) {
				tempAngleXZ = (float) (Math.atan2(elm.y[0][i + 1] * scale
						- elm.y[0][i] * scale, elm.x[i + 1] * elm.dispLScale
						- elm.x[i] * elm.dispLScale));
				tempAngleXY = (float) (Math.atan2(elm.y[1][i + 1] * scale
						- elm.y[1][i] * scale, elm.x[i + 1] * elm.dispLScale
						- elm.x[i] * elm.dispLScale));

			}
			p.vertex(
					elm.x[i] * elm.dispLScale - beamHeight
							* PApplet.sin(tempAngleXZ) + beamWidth
							* PApplet.sin(tempAngleXY), elm.y[1][i] * scale
							- beamWidth * PApplet.cos(tempAngleXY), elm.y[0][i]
							* scale + beamHeight * PApplet.cos(tempAngleXZ));
			p.vertex(
					elm.x[i] * elm.dispLScale - beamHeight
							* PApplet.sin(tempAngleXZ) - beamWidth
							* PApplet.sin(tempAngleXY), elm.y[1][i] * scale
							+ beamWidth * PApplet.cos(tempAngleXY), elm.y[0][i]
							* scale + beamHeight * PApplet.cos(tempAngleXZ));

		}
		p.endShape(PConstants.CLOSE);
	}

	public float getIxz() {
		return Ixz;
	}

	public void setIxz(float ixz) {
		Ixz = ixz;
	}

	public float getIxy() {
		return Ixy;
	}

	public void setIxy(float ixy) {
		Ixy = ixy;
	}

	public float getJ() {
		return J;
	}

	public void setJ(float j) {
		J = j;
	}

	public float getA() {
		return A;
	}

	public void setA(float a) {
		A = a;
	}

	public float getBeamHeight() {
		return beamHeight;
	}

	public void setBeamHeight(float beamHeight) {
		this.beamHeight = beamHeight;
	}

	public float getBeamWidth() {
		return beamWidth;
	}

	public void setBeamWidth(float beamWidth) {
		this.beamWidth = beamWidth;
	}

}
