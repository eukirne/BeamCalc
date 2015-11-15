/**
 * (./) Topology_Examples.java v0.1 05/09/2011
 * @author  Enrique Ramos Melgar. 
 * http://www.esc-studio.com
 * 
 * This is an example application to test the BeamCalc library.
 * http://www.processing.org
 */

import processing.core.PApplet;
import processing.core.PVector;
import BeamCalc.Structure;

@SuppressWarnings("serial")
public class Topology_Examples extends PApplet {
	// Add main function to create Application
	public static void main(String args[]) {
		PApplet.main(new String[] { "Topology_Examples" });
	}

	// Structure class
	Structure structure;

	public void setup() {

		size(800, 700, OPENGL);
		smooth();

		// Initialise Structure
		// The second parameter defines the number of points
		// for internal stress analysis
		structure = new Structure(this, 7);
		initialise(0);

	}

	public void initialise(int topology) {
		// Initialise Topology
		switch (topology) {
		case 0:
			arup00();
			break;
		case 1:
			arup01();
			break;
		case 2:
			arup02();
			break;
		case 3:
			frame02();
			break;
		case 4:
			frame03();
			break;
		case 5:
			frame04();
			break;
		}
		structure.render(false);
	}

	public void draw() {

		// This conditional takes values from the
		// structure class, and initialises the different
		// topologies when the button is pressed
		if (structure.changeTopology) {
			structure.reset();
			initialise(structure.topology);
			structure.changeTopology = false;
		}

		background(255);

		// Evaluate Structure
		structure.eval();

		// Draw Structure
		structure.draw(this);

	}

	private void arup00() {

		// Three element structure, pinned nodes
		PVector pt1 = new PVector(0, 0, 0);
		PVector pt2 = new PVector(0, 0, 6000);
		PVector pt3 = new PVector(6000, 0, 6000);
		PVector pt4 = new PVector(9000, 0, 0);
		structure.addBeam(pt1, pt2, 1, 0, 0);
		structure.addBeam(pt2, pt3, 0, 0, 0);
		structure.addBeam(pt3, pt4, 0, 1, 0);

	}

	private void arup01() {

		// Three element structure, clamped nodes
		PVector pt1 = new PVector(0, 0, 0);
		PVector pt2 = new PVector(0, 0, 5000);
		PVector pt3 = new PVector(5000, 0, 5000);
		PVector pt4 = new PVector(5000, 0, 0);

		structure.addBeam(pt1, pt2, 3, 0, 0);
		structure.addBeam(pt2, pt3, 0, 0, 0);
		structure.addBeam(pt3, pt4, 0, 3, 0);

	}

	private void arup02() {

		// Complex two-dimensional structure
		PVector pt1 = new PVector(0, 0, 0);
		PVector pt2 = new PVector(0, 0, 6000);
		PVector pt3 = new PVector(0, 0, 12000);
		PVector pt4 = new PVector(0, 0, 18000);
		PVector pt5 = new PVector(6000, 0, 0);
		PVector pt6 = new PVector(6000, 0, 6000);
		PVector pt7 = new PVector(6000, 0, 12000);
		PVector pt8 = new PVector(6000, 0, 18000);
		PVector pt9 = new PVector(12000, 0, 0);
		PVector pt10 = new PVector(12000, 0, 6000);
		PVector pt11 = new PVector(12000, 0, 12000);
		PVector pt12 = new PVector(12000, 0, 18000);
		PVector pt13 = new PVector(18000, 0, 0);
		PVector pt14 = new PVector(18000, 0, 6000);
		PVector pt15 = new PVector(18000, 0, 12000);
		PVector pt16 = new PVector(18000, 0, 18000);
		PVector pt17 = new PVector(24000, 0, 6000);
		PVector pt18 = new PVector(24000, 0, 12000);
		PVector pt19 = new PVector(24000, 0, 18000);

		structure.addBeam(pt1, pt2, 3, 0, 0);
		structure.addBeam(pt2, pt3, 0, 0, 0);
		structure.addBeam(pt3, pt4, 0, 0, 0);
		structure.addBeam(pt5, pt6, 3, 0, 0);
		structure.addBeam(pt6, pt7, 0, 0, 0);
		structure.addBeam(pt7, pt8, 0, 0, 0);
		structure.addBeam(pt9, pt10, 3, 0, 0);
		structure.addBeam(pt10, pt11, 0, 0, 0);
		structure.addBeam(pt11, pt12, 0, 0, 0);
		structure.addBeam(pt13, pt14, 1, 0, 0);
		structure.addBeam(pt14, pt15, 0, 0, 0);
		structure.addBeam(pt15, pt16, 0, 0, 0);
		structure.addBeam(pt18, pt19, 0, 0, 0);
		structure.addBeam(pt2, pt6, 0, 0, 0);
		structure.addBeam(pt3, pt7, 0, 0, 0);
		structure.addBeam(pt4, pt8, 0, 0, 0);
		structure.addBeam(pt6, pt10, 0, 0, 0);
		structure.addBeam(pt7, pt11, 0, 0, 0);
		structure.addBeam(pt8, pt12, 0, 0, 0);
		structure.addBeam(pt10, pt14, 0, 0, 0);
		structure.addBeam(pt11, pt15, 0, 0, 0);
		structure.addBeam(pt12, pt16, 0, 0, 0);
		structure.addBeam(pt14, pt17, 0, 0, 0);
		structure.addBeam(pt15, pt18, 0, 0, 0);
		structure.addBeam(pt16, pt19, 0, 0, 0);

	}

	private void frame02() {

		// Three-dimensional Structure
		PVector pt1 = new PVector(0, 0, 0);
		PVector pt2 = new PVector(0, 0, 5000);
		PVector pt3 = new PVector(5000, 0, 5000);
		PVector pt4 = new PVector(5000, 0, 0);
		PVector pt5 = new PVector(0, 5000, 5000);
		PVector pt6 = new PVector(0, 5000, 0);
		PVector pt7 = new PVector(5000, 5000, 5000);
		PVector pt8 = new PVector(5000, 5000, 0);
		PVector pt9 = new PVector(-5000, 5000, 5000);
		PVector pt10 = new PVector(0, 0, 10000);

		structure.addBeam(pt1, pt2, 3, 0, 0);
		structure.addBeam(pt2, pt3, 0, 0, 0);
		structure.addBeam(pt3, pt4, 0, 3, 0);
		structure.addBeam(pt5, pt2, 0, 0, 1);
		structure.addBeam(pt5, pt6, 0, 3, 1);
		structure.addBeam(pt3, pt7, 0, 0, 1);
		structure.addBeam(pt7, pt8, 0, 3, 1);
		structure.addBeam(pt5, pt7, 0, 0, 0);
		structure.addBeam(pt9, pt5, 0, 0, 2);
		structure.addBeam(pt2, pt10, 0, 0, 0);

	}

	private void frame03() {

		// Three-dimensional irregular frame
		PVector[][][] myVec;

		int dim = 4;
		myVec = new PVector[dim][dim][dim];
		float size = 5000;
		float noise = 800;
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				for (int k = 0; k < dim; k++) {
					myVec[i][j][k] = new PVector(size * i
							+ (float) Math.random() * noise, size * j
							+ (float) Math.random() * noise, size * k
							+ (float) Math.random() * noise);
				}
			}
		}

		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				for (int k = 0; k < dim; k++) {

					if (k != 0) {
						if (i < dim - 1)
							structure.addBeam(myVec[i][j][k],
									myVec[i + 1][j][k], 0, 0, 0);
						if (j < dim - 1)
							structure.addBeam(myVec[i][j][k],
									myVec[i][j + 1][k], 0, 0, 2);
					}
					if (k < dim - 1)
						structure.addBeam(myVec[i][j][k], myVec[i][j][k + 1],
								0, 0, 1);

				}
			}
		}

		for (int i = 0; i < structure.nodes.size(); i++) {
			if (structure.nodes.get(i).pos.z < 1000) {
				structure.nodes.get(i).linked = true;
				structure.nodes.get(i).clamped = true;
			}
		}
	}

	private void frame04() {

		// Elliptical Paraboloid
		PVector[][] vec;

		int dim = 5;
		vec = new PVector[2 * dim][2 * dim];
		float size = 2000;
		float noise = 0;
		float height = 100;
		for (int i = 0; i < 2 * dim; i++) {
			for (int j = 0; j < 2 * dim; j++) {
				vec[i][j] = new PVector(size * (i - dim + 1)
						+ (float) Math.random() * noise, size * (j - dim + 1)
						+ (float) Math.random() * noise, height * (i - dim + 1)
						* (i - dim + 1) - height * (j - dim + 1)
						* (j - dim + 1) + (float) Math.random() * noise);
			}
		}

		for (int i = 0; i < 2 * dim; i++) {
			for (int j = 0; j < 2 * dim; j++) {
				if (i < 2 * dim - 1 && j != 0 && j != 2 * dim - 1)
					structure.addBeam(vec[i][j], vec[i + 1][j], 0, 0, 0);

				if (j < 2 * dim - 1)
					if (j == 0) {
						structure.addBeam(vec[i][j], vec[i][j + 1], 1, 0, 1);
					} else if (j == 2 * dim - 2) {
						structure.addBeam(vec[i][j], vec[i][j + 1], 0, 1, 1);
					} else {
						structure.addBeam(vec[i][j], vec[i][j + 1], 0, 0, 1);
					}
			}
		}

	}

}
