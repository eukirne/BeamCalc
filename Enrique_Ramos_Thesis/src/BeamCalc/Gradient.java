/**
 * (./) Gradient.java v0.1 05/09/2011
 * @author Enrique Ramos Melgar
 * http://www.esc-studio.com
 *
 * THIS LIBRARY IS RELEASED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE
 * http://creativecommons.org/licenses/by/3.0/
 * http://www.processing.org/
 * 
 */

package BeamCalc;

/**
 * Gradient Class, used to create gradients containing the 
 * whole colour spectrum
 * @author Enrique Ramos
 *
 */
public class Gradient {

	private float min;
	private float max;

	public Gradient() {

	}
	
	public Gradient(float min, float max){
		this.min = min;
		this.max = max;
	}

	public int[] getColor(float t) {
		int[] col = new int[3];
		//float t=min + (param)*(max - min);

		if (t < min + (max - min) * 0.25) {
			col[0] = 0;
			col[1] = (int) ((255 * 4 / (max - min) * (t - min)));
			col[2] = 255;
		}
		if (t >= min + (max - min) * 0.25 && t < min + (max - min) * 0.5) {
			col[0] = 0;
			col[1] = 255;
			col[2] = (int) (255 - (255 * 4 / (max - min) * (t
					- min - (max - min) * 0.25)));
		}
		if (t >= min + (max - min) * 0.5 && t < min + (max - min) * 0.75) {
			col[0] = (int) ((255 * 4 / (max - min) * (t - min - (max - min) * 0.5)));
			col[1] = 255;
			col[2] = 0;
		}
		if (t >= min + (max - min) * 0.75) {
			col[0] = 255;
			col[1] = (int) (255 - (255 * 4 / (max - min) * (t
					- min - (max - min) * 0.75)));
			col[2] = 0;
		}

		return col;
	}

	public void setBounds(float min, float max){
		this.min = min;
		this.max = max;
	}
}
