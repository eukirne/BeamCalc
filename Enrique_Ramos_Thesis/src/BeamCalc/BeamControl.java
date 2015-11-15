/**
 * (./) BeamControl.java v0.1 05/09/2011
 * @author Enrique Ramos Melgar
 * http://www.esc-studio.com
 *
 * THIS LIBRARY IS RELEASED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE
 * http://creativecommons.org/licenses/by/3.0/
 * http://www.processing.org/
 */
package BeamCalc;

import java.text.DecimalFormat;

import javax.swing.JFrame;

import processing.core.PApplet;
import processing.core.PFont;

@SuppressWarnings("serial")
/**
 * This class creates a new instance of JFrame containing a 
 * Processing application. The element's material and geometric 
 * properties can be visualised and modified here.
 * @author Enrique Ramos
 */
public class BeamControl extends JFrame {

	Structure world;
	DisplayData s;

	public BeamControl(Structure world) {
		this.world = world;
		this.setTitle("Section Properties");
		setBounds(world.p.frame.getLocation().x + world.p.frame.getWidth() + 5,
				world.p.frame.getLocation().y, 200,
				world.p.frame.getHeight() / 2);
		s = new DisplayData(world, this);
		this.add(s);
		s.init();
		// this.setResizable(false);
		// this.setFocusable(false);
		// setUndecorated(true);
		setVisible(true); // was show();

	}

	public void draw() {
		s.draw();
	}

	public class DisplayData extends PApplet {

		// Orbit class
		Structure world;
		BeamControl p;
		int frameNum = 0;
		boolean action = false;

		PFont fontA;
		int textY;

		int mouseDrag;
		int previousDrag;

		DecimalFormat df = new DecimalFormat("#0.000");

		public boolean display = false;
		private float sectionWidth;
		private float sectionHeight;
		private int mouse_x;
		private int mouse_y;
		private float prevSectionWidth;
		private float prevSectionHeight;
		private float sectionIxy;
		private float sectionIxz;
		private float prevSectionIxy;
		private float prevSectionIxz;

		DisplayData(Structure world, BeamControl p) {
			this.world = world;
			this.p = p;
		}

		public void setup() {
			size(p.getWidth(), p.getHeight());
			smooth();
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

				sectionWidth = world.focusElms.get(0).section.getBeamWidth();
				sectionHeight = world.focusElms.get(0).section.getBeamHeight();
				sectionIxy = world.focusElms.get(0).section.getIxy();
				sectionIxz = world.focusElms.get(0).section.getIxz();
				
				textAlign(LEFT);
				fill(100, 0, 0);
				translate(width / 8, 0);
					translate(0, (0 + 1) * height / 4);
					textY += (0 + 1) * height / 4;
					fill(100, 0, 0);
					textAlign(LEFT);
					text("ELEMENT " + world.focusElms.get(0).id, -20, -55);
					
					text("Element Section ", -20, -35);
					text("Ixy = "+ sectionIxy, -20, -15);
					text("Ixz = "+ sectionIxz, -20, -5);

					stroke(0);
					fill(200);
					rect(0, 0, sectionWidth, sectionHeight);


			} catch (Exception exc) {

			}
		}

		public void mousePressed() {
			mouse_x = mouseX;
			mouse_y = mouseY;
			prevSectionWidth = sectionWidth;
			prevSectionHeight = sectionWidth;
			prevSectionIxy = sectionIxy;
			prevSectionIxz = sectionIxz;

		}

		public void mouseClicked() {

		}

		public void mouseDragged() {
			if (mouseButton == PApplet.LEFT) {
				world.focusElms.get(0).section.setBeamWidth(prevSectionWidth + mouseX - mouse_x);
				world.focusElms.get(0).section.setBeamHeight(prevSectionHeight + mouseY - mouse_y);
				world.focusElms.get(0).section.setIxy(prevSectionIxy + 100000*(mouseX - mouse_x));
				world.focusElms.get(0).section.setIxz(prevSectionIxz + 100000*(mouseY - mouse_y));

			}if (mouseButton == PApplet.RIGHT) {

			}
		}

	}
}
