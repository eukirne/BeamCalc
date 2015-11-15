/**
 * (./) orbit.java v0.1 28/05/2011
 * (by) Enrique Ramos Melgar
 * http://www.esc-studio.com
 * 
 * This is a basic orbit library for Procesing based on 
 * code from Przemek Jaworsky - http://www.jawordesign.com
 *
 * http://www.processing.org/
 *
 * THIS LIBRARY IS RELEASED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE
 * http://creativecommons.org/licenses/by/3.0/
 */

package orbit.library;

import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.io.File;
import java.io.PrintWriter;
import processing.core.*;

/**
 * This Class creates an orbit object that allows the use of the 
 * muose to control the three-dimensional rotation of the space.
 * @author Enrique Ramos
 *
 */
public class Orbit {

	private PApplet parent;
	// Camera Variables
	String[] viewPoint = new String[7];
	String viewPointFile;

	/*
	 * Orbit Parameters as stored in the csv file orbitParams[0] = alpha1
	 * orbitParams[1] = alpha2 orbitParams[2] = distance orbitParams[3] =
	 * translate1 orbitParams[4] = translate2 orbitParams[5] = translate3
	 * orbitParams[6] = fov
	 */
	private float[] orbitParams = new float[7];

	private float panScale = 50;
	private float CSScale = 100;

	private float alpha1_, alpha2_;
	private int mouse_x, mouse_y;
	private boolean retPers = true;
	private int numRot = 0;
	public boolean drawCS = true;
	private PrintWriter output;
	private File f;
	boolean milimeters;

	private boolean shiftMode;

	/**
	 * Main constructor. It takes a Processing Application and
	 * a number for the file containing the saved views
	 * @param p
	 * @param numView
	 */
	public Orbit(PApplet p, int numView) {
		this.parent = p;
		parent.registerDispose(this);
		parent.registerMouseEvent(this);
		parent.registerKeyEvent(this);

		mouseWheel();

		// Creates a file containing the saved view settings.
		viewPointFile = "orbitSet_" + numView + ".csv";

		// Check if there is a saved view in the data folder of the sketch, and
		// create a file if negative
		f = new File(p.dataPath(viewPointFile));
		if (!f.exists()) {
			output = parent.createWriter(p.dataPath(viewPointFile));

			if (milimeters) {
				output.println(-6.81999);
				output.println(7.829951);
				output.println(7.333359);
				output.println(-24.448734);
				output.println(11.007046);
				output.println(-948.3287);
				output.println(0.7843585);
			} else {
				output.println(0.15);
				output.println(1.15);
				output.println(0.5);
				output.println(-9.0);
				output.println(4.5);
				output.println(15.5);
				output.println(0.8);
			}
			output.flush(); // Write the remaining data
			output.close(); // Finish the file
		}
	}

	// -----------------------------------------------------------------------
	// PushOrbit and PopOrbit methods

	public float getPanScale() {
		return panScale;
	}

	public void setPanScale(float panScale) {
		this.panScale = panScale;
	}

	public void pushOrbit(PApplet win) {
		pushOrbit(win.g);
	}
	
	public void popOrbit(PApplet win) {
		popOrbit(win.g);
	}
	
	public void pushOrbit(PGraphics window) {

		window.lights();
		window.camera((float) (parent.width / 2.0),
				(float) (parent.height / 2.0), (float) (orbitParams[2]
						* (parent.height / 2.0) / Math
						.tan(Math.PI * 60.0 / 360.0)),
				(float) (parent.width / 2.0), (float) (parent.height / 2.0),
				(float) (-1 * (parent.height / 2.0) / Math
						.tan(Math.PI * 60.0 / 360.0)), 0f, 1f, 0f);
		window.perspective(
				orbitParams[6],
				(((float) parent.width) / ((float) parent.height)),
				(float) ((orbitParams[2]
						* ((parent.height / 2.0) / Math
								.tan(Math.PI * 60.0 / 360.0)) / 10.0f)),
				(float) (((parent.height / 2.0) / Math
						.tan(Math.PI * 60.0 / 360.0)) * 500.0));

		alpha2_ = (float) (orbitParams[1] + Math.sin(parent.frameCount / 130.0) * 0.00008);
		alpha1_ = (float) (orbitParams[0] + Math.sin(parent.frameCount / 123.0) * 0.00008);

		window.translate(parent.width / 2, parent.height / 2);
		window.rotateX(alpha2_);
		window.rotateZ(alpha1_);
		window.translate(orbitParams[3], orbitParams[4], orbitParams[5]);
		window.scale(1, -1, 1);

		if (retPers) {
			retrievePerspectiveSettings();
		}

		if (drawCS)
			drawCS(window,CSScale);
	}

	public void popOrbit(PGraphics window) {
		window.camera();
		window.perspective();
	}

	public float getCSScale() {
		return CSScale;
	}

	public void setCSScale(float cSScale) {
		CSScale = cSScale;
	}

	void drawCS(PGraphics window, float cSScale2) {
		// Draw Coordinate System
		window.strokeWeight(1);
		window.stroke(255, 0, 0);
		window.line(0, 0, 0, cSScale2, 0, 0);
		window.stroke(0, 255, 0);
		window.line(0, 0, 0, 0, cSScale2, 0);
		window.stroke(0, 0, 255);
		window.line(0, 0, 0, 0, 0, cSScale2);
		window.stroke(0);
	}

	// -----------------------------------------------------------------------
	// File Operations

	private void savePerspectiveSettings() {

		for (int i = 0; i < orbitParams.length; i++) {

			viewPoint[i] = String.valueOf(orbitParams[i]);
		}

		PApplet.saveStrings(f, viewPoint);
	}

	private void retrievePerspectiveSettings() {

		String lines[] = PApplet.loadStrings(f);

		if (f.exists()) {

			// Set the number of frame for perspective change animation
			if (numRot < 50) {

				for (int i = 0; i < orbitParams.length; i++) {

					orbitParams[i] = (Float.valueOf(lines[i]).floatValue() + orbitParams[i]) * 0.5f;
					numRot++;
				}
			} else {
				for (int i = 0; i < orbitParams.length; i++) {
					orbitParams[i] = Float.valueOf(lines[i]).floatValue();
				}

				retPers = false;
				numRot = 0;
			}
		}
	}

	// -----------------------------------------------------------------------
	// Keys and Mouse

	public void keyEvent(KeyEvent e) {

		if (parent.key == 'p') {
			savePerspectiveSettings();
		}

		if (parent.key == 'o') {
			retPers = true;
		}
		
		if (e.isShiftDown()) {
			shiftMode = true;
		} else {
			shiftMode = false;
		}

	}

	public void mouseEvent(MouseEvent event) {
		int x = event.getX();
		int y = event.getY();

		switch (event.getID()) {
		case MouseEvent.MOUSE_PRESSED:
			mouse_x = x;
			mouse_y = y;
			break;
		case MouseEvent.MOUSE_RELEASED:
			// do something for mouse released
			break;
		case MouseEvent.MOUSE_CLICKED:
			// do something for mouse clicked
			break;
		case MouseEvent.MOUSE_DRAGGED:
			if (!shiftMode) {
				if (parent.mouseButton == PConstants.LEFT) {
					orbitParams[0] += (mouse_x - parent.mouseX) / 100.0;
					orbitParams[1] += (mouse_y - parent.mouseY) / 100.0;
				}

				if (parent.mouseButton == PConstants.RIGHT) {
					float msx = (mouse_x - parent.mouseX);
					float msy = (mouse_y - parent.mouseY)
							* PApplet.cos(orbitParams[1]);
					float msz = (mouse_y - parent.mouseY);
					orbitParams[3] -= (Math.cos(orbitParams[0]) * msx + Math
							.sin(orbitParams[0]) * msy)
							* panScale;
					orbitParams[4] -= (-Math.sin(orbitParams[0]) * msx + Math
							.cos(orbitParams[0]) * msy)
							* panScale;
					orbitParams[5] += (Math.sin(orbitParams[1]) * msz)
							* panScale;
				}

				if (parent.mouseButton == PConstants.CENTER) {
					orbitParams[6] = orbitParams[6]
							* (parent.mouseX - mouse_x + parent.width)
							/ parent.width;
					float msy = (-mouse_y + parent.mouseY);
					orbitParams[2] = orbitParams[2] * (msy + 20) / 20;
				}

				mouse_x = parent.mouseX;
				mouse_y = parent.mouseY;
			}

			break;

		case MouseEvent.MOUSE_MOVED:
			break;
		}
	}

	private void mouseWheel() {
		MouseWheelListener wheelListener = new MouseWheelListener() {
			@Override
			public void mouseWheelMoved(MouseWheelEvent e) {
				mouseWheel(e.getWheelRotation());
			}
		};
		parent.addMouseWheelListener(wheelListener);
	}

	private void mouseWheel(int delta) {

		orbitParams[2] = orbitParams[2] * (delta + 20) / 20;

	}

	public void dispose() {
		// anything in here will be called automatically when
		// the parent applet shuts down. for instance, this might
		// shut down a thread used by this library.
	}

}
