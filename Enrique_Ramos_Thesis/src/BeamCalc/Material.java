/**
 * (./) Material.java v0.1 05/09/2011
 * @author Enrique Ramos Melgar
 * http://www.esc-studio.com
 *
 * THIS LIBRARY IS RELEASED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE
 * http://creativecommons.org/licenses/by/3.0/
 * http://www.processing.org/
 */

package BeamCalc;

/**
 * Class containing the material properties of the element
 * @author Enrique Ramos
 *
 */
public class Material {

	// Material Properties
	private float E; // modulus of elasticity (N/mm2)
	private float G; // Shear Modulus (N/mm2)
	private float D; // Density (T/mm2)
	
	public Material(){
		
		this.E = 205000; // N/mm2
		this.G = 78850; // N/mm2
		this.D = 0.0000785f; // Density (N/mm3)
	}

	public float getE() {
		return E;
	}

	public void setE(float e) {
		E = e;
	}

	public float getG() {
		return G;
	}

	public void setG(float g) {
		G = g;
	}

	public float getD() {
		return D;
	}

	public void setD(float d) {
		D = d;
	}
	
}
