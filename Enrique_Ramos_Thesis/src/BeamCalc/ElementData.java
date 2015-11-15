/**
 * (./) ElementData.java v0.1 05/09/2011
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
import processing.core.PConstants;
import processing.core.PFont;

@SuppressWarnings("serial")
/**
 * This class creates a java JFrame containing a Processing application
 * in which the data of the selected elements is displayed, and can be modified 
 * by the user
 * @author Enrique Ramos
 */
public class ElementData extends JFrame {

	Structure world;
	DisplayData s;

	public ElementData(Structure world) {
		this.world = world;
		this.setTitle("Internal Stresses");
		setBounds(world.p.frame.getLocation().x - 255,
				world.p.frame.getLocation().y+world.p.frame.getHeight()/2+5
				, 250, world.p.frame.getHeight()/2);
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
		ElementData p;
		int frameNum = 0;
		boolean action = false;

		PFont fontA;
		int textY;

		int mouseDrag;
		int previousDrag;

		DecimalFormat df = new DecimalFormat("#0.000");

		public boolean display = false;

		DisplayData(Structure world, ElementData p) {
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
				
				textAlign(LEFT);
				fill(100, 0, 0);		
				translate(width/8,0);
				for (int i = 0; i < world.focusElms.size(); i++) {
									
					pushMatrix();
					translate(0, (4*i+1)*height / 4);
					textY+=(4*i+1)*height / 4;
					fill(100, 0, 0);
					textAlign(LEFT);
					text("ELEMENT " + world.focusElms.get(i).id, -20,-55);
					drawMoment(world.focusElms.get(i), 1);
					translate(0, height / 4);
					drawShear(world.focusElms.get(i), 1);
					translate(0, height / 4);
					drawDisp(world.focusElms.get(i), 10);
					popMatrix();
				}
				
			} catch (Exception exc) {

			}
		}

		private void drawDisp(Element fElm, float scale) {

			float yScale = -0.05f * height / fElm.getMaxDisp();
			float xScale = 0.8f * width / fElm.L;
			
			noFill();

			stroke(100);
			line(0,0,fElm.L*xScale,0);
			strokeWeight(2);
			stroke(0);
			for (int i = 0; i < fElm.numPoints - 1; i++) {
				line(fElm.x[i] * xScale, fElm.y[1][i] * yScale,
						fElm.x[i + 1] * xScale, fElm.y[1][i + 1]
								* yScale);
			}
			
			translate(0, 50);

			strokeWeight(1);
			line(0,0,fElm.L*xScale,0);
			strokeWeight(2);
			stroke(0);
			for (int i = 0; i < fElm.numPoints - 1; i++) {
				line(fElm.x[i] * xScale, fElm.y[0][i] * yScale,
						fElm.x[i + 1] * xScale, fElm.y[0][i + 1] * yScale);
			}

			strokeWeight(1);

		}

		private void drawMoment(Element fElm, float scale) {

			float yScale = -0.1f * height / fElm.getMaxMoment();
			float xScale = 0.8f * width / fElm.L;

			fill(0);
			textAlign(LEFT);
			
			text("MOMENTS",0,-0.1f*height);
			text("Mxy",-25,0);
			text(fElm.aMoment.z/1000,0,0);
			textAlign(RIGHT);
			text(fElm.bMoment.z/1000,fElm.L*xScale,0);
			
			stroke(150, 50, 0);
			fill(100, 50);
			beginShape(PConstants.QUAD_STRIP);
			for (int i = 0; i < fElm.numPoints; i++) {
				vertex(fElm.x[i] * xScale, 0);
				vertex(fElm.x[i] * xScale, fElm.m[1][i] * yScale);
			}
			endShape(PConstants.CLOSE);

			translate(0, 50);
			
			fill(0);
			textAlign(LEFT);
			text("Mxz",-25,0);
			text(fElm.aMoment.y/1000,0,0);
			textAlign(RIGHT);
			text(fElm.bMoment.y/1000,fElm.L*xScale,0);
			
			
			stroke(150, 50, 0);
			// p.noStroke();
			fill(100, 50);
			// p.noFill();
			beginShape(PConstants.QUAD_STRIP);
			for (int i = 0; i < fElm.numPoints; i++) {
				vertex(fElm.x[i] * xScale, 0);
				vertex(fElm.x[i] * xScale, fElm.m[0][i] * yScale);
			}
			endShape(PConstants.CLOSE);
		}
		
		private void drawShear(Element fElm, float scale) {

			float yScale = -0.1f * height / fElm.getMaxShear();
			float xScale = 0.8f * width / fElm.L;

			fill(0);
			textAlign(LEFT);
			text("SHEAR",0,-0.1f*height);
			text("Vxy",-25,0);
			text(fElm.aReact.z/1000,0,0);
			textAlign(RIGHT);
			text(fElm.bReact.z/1000,fElm.L*xScale,0);
			
			stroke(50, 100, 0);
			fill(100, 50);
			beginShape(PConstants.QUAD_STRIP);
			for (int i = 0; i < fElm.numPoints; i++) {
				vertex(fElm.x[i] * xScale, 0);
				vertex(fElm.x[i] * xScale, fElm.v[1][i] * yScale);
			}
			endShape(PConstants.CLOSE);

			translate(0, 50);
			
			fill(0);
			textAlign(LEFT);
			text("Vxz",-25,0);
			text(fElm.aReact.y/1000,0,0);
			textAlign(RIGHT);
			text(fElm.bReact.y/1000,fElm.L*xScale,0);
			
			
			stroke(50, 100, 0);
			// p.noStroke();
			fill(100, 50);
			// p.noFill();
			beginShape(PConstants.QUAD_STRIP);
			for (int i = 0; i < fElm.numPoints; i++) {
				vertex(fElm.x[i] * xScale, 0);
				vertex(fElm.x[i] * xScale, fElm.v[0][i] * yScale);
			}
			endShape(PConstants.CLOSE);
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
