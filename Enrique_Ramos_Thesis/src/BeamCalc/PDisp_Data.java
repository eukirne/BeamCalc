/**
 * (./) PDisp_Data.java v0.1 05/09/2011
 * @author Enrique Ramos Melgar
 * http://www.esc-studio.com
 *
 * THIS LIBRARY IS RELEASED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE
 * http://creativecommons.org/licenses/by/3.0/
 * http://www.processing.org/
 * 
 */

package BeamCalc;

import java.text.DecimalFormat;

import javax.swing.JFrame;

import processing.core.PApplet;
import processing.core.PFont;

@SuppressWarnings("serial")
/**
 * This class creates a java JFrame containing a Processing application
 * in which the data of the selected elements is displayed.
 * @author Enrique Ramos
 */
public class PDisp_Data extends JFrame {

	Structure world;
	DisplayData s;

	public PDisp_Data(Structure world) {
		this.world = world;
		this.setTitle("Structure Data");
		setBounds(world.p.frame.getLocation().x - 255,
				world.p.frame.getLocation().y, 250, world.p.frame.getHeight()/2);
		s = new DisplayData(world);
		this.add(s);
		s.init();
//		this.setResizable(false);
//		this.setFocusable(false);
//		setUndecorated(true);
		setVisible(true); // was show();

	}

	public void draw() {
		s.draw();
	}

	public class DisplayData extends PApplet {

		// Orbit class
		Structure world;
		int frameNum = 0;
		boolean action = false;

		PFont fontA;
		int textY;

		int mouseDrag;
		int previousDrag;

		DecimalFormat df = new DecimalFormat("#0.000");

		public boolean display = false;

		DisplayData(Structure world) {
			this.world = world;
		}

		public void setup() {
			size(250, world.p.frame.getHeight()/2-5);
			// noLoop();
			fontA = createFont("ArialMT-10.vlw", 10);
			textFont(fontA, 12);
		}

		public void draw() {

			background(240);
			stroke(0);
			fill(0);
			textFont(fontA, 10);

			mouseDrag = constrain(mouseDrag, -abs(height - textY), 0);
			translate(0, mouseDrag);

			textY = 20;

			try {
				fill(100, 0, 0);
				text("Structure Data ", 10, textY);
				textY += 15;
				fill(0);
				text("Number of Elements: " + world.elms.size(), 10, textY);
				textY += 15;
				text("Number of Nodes: " + world.nodes.size(), 10, textY);
				textY += 30;
				
				int end = millis();
				int timeSpent = (end) / 1000; // Seconds spent
				int h = timeSpent / 3600;
				timeSpent -= h * 3600;
				int m = timeSpent / 60;
				timeSpent -= m * 60;
				
				text("Calculation Time: " +	h + " h, " + m + " m, " + timeSpent + " s", 10, textY);
				textY += 15;
				text("Iterations: " + world.iterations, 10, textY);
				textY += 30;

				for (int i = 0; i < world.focusElms.size(); i++) {
					fill(100, 0, 0);
					text("Element " + world.focusElms.get(i).id, 10, textY);
					textY += 15;
					fill(0);
					text("	Element Angle X:  "
							+ df.format(Math.toDegrees(world.focusElms.get(i).aAngles.x)),
							10, textY);
					textY += 15;
					text("	Element Angle Y:  "
							+ df.format(Math.toDegrees(world.focusElms.get(i).aAngles.y)),
							10, textY);
					textY += 15;
					text("	Element Angle Z:  "
							+ df.format(Math.toDegrees(world.focusElms.get(i).aAngles.z)),
							10, textY);
					textY += 15;
					text("	Element Initial Length:  "
							+ df.format(world.focusElms.get(i).initialL), 10,
							textY);
					textY += 15;
					text("	Element Length:  "
							+ df.format(world.focusElms.get(i).L), 10, textY);
					textY += 15;
					text("	Node A cVel:  " + world.focusElms.get(i).a.damping[0].mag(), 10,
							textY);
					textY += 15;
					text("	Node A PosStep:  " + world.focusElms.get(i).a.posStep, 10,
							textY);
					textY += 15;
					text("	Node B id:  " + world.focusElms.get(i).b.id, 10,
							textY);
					textY += 15;
					text("	Node B PosStep:  " + world.focusElms.get(i).b.posStep, 10,
							textY);
					textY += 15;
					text("	Element axial:  "
							+ df.format(world.focusElms.get(i).axial), 10,
							textY);
					textY += 15;
					text("	Element WeightReact:  "
							+ df.format(world.focusElms.get(i).weightAxial),
							10, textY);
					textY += 15;
					text("	Node A Rotation:  "
							+ df.format(Math.toDegrees(world.focusElms.get(i).a.rot.x)),
							10, textY);
					textY += 15;
					text("	Node B Rotation:  "
							+ df.format(Math.toDegrees(world.focusElms.get(i).b.rot.x)),
							10, textY);
					textY += 15;
					text("	Node A X Displacement:  "
							+ df.format(world.focusElms.get(i).a.dispVec.x),
							10, textY);
					textY += 15;
					text("	Node A Y Displacement:  "
							+ df.format(world.focusElms.get(i).a.dispVec.y),
							10, textY);
					textY += 15;
					text("	Node A Z Displacement:  "
							+ df.format(world.focusElms.get(i).a.dispVec.z),
							10, textY);
					textY += 15;
					text("	Node B X Displacement:  "
							+ df.format(world.focusElms.get(i).b.dispVec.x),
							10, textY);
					textY += 15;
					text("	Node B Y Displacement:  "
							+ df.format(world.focusElms.get(i).b.dispVec.y),
							10, textY);
					textY += 15;
					text("	Node B Z Displacement:  "
							+ df.format(world.focusElms.get(i).b.dispVec.z),
							10, textY);
					textY += 15;
					text("	Moment A:  "
							+ df.format(world.focusElms.get(i).mm1[0] / 1000)
							+ "  mm1: "
							+ df.format(world.focusElms.get(i).m[0][0]), 10,
							textY);
					textY += 15;
					text("	Moment B:  "
							+ df.format(world.focusElms.get(i).mm2[0] / 1000)
							+ "  mm2: "
							+ df.format(world.focusElms.get(i).m[0][world.focusElms
									.get(i).numPoints - 1]), 10, textY);
					textY += 15;
					text("	Moment A node:  "
							+ df.format(world.focusElms.get(i).a.momentVec.x),
							10, textY);
					textY += 15;
					text("	Moment B node:  "
							+ df.format(world.focusElms.get(i).b.momentVec.x),
							10, textY);
					textY += 15;
					text("	Reaction A node:  "
							+ df.format(world.focusElms.get(i).rr1[0]), 10,
							textY);
					textY += 15;
					text("	Reaction B node:  "
							+ df.format(world.focusElms.get(i).rr2[0]), 10,
							textY);
					textY += 30;
				}
			} catch (Exception exc) {

			}
		}

		public void mousePressed() {
			previousDrag = mouseY - mouseDrag;
		}

		public void mouseReleased() {
			previousDrag = mouseDrag;
		}

		public void mouseDragged() {
			mouseDrag = mouseY - previousDrag;
		}
	}
}
