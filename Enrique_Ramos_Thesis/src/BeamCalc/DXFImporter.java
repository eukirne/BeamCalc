/**
 * (./) DXFImporter.java v0.1 05/09/2011
 * @author Enrique Ramos Melgar
 * http://www.esc-studio.com
 *
 * THIS LIBRARY IS RELEASED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE
 * http://creativecommons.org/licenses/by/3.0/
 * http://www.processing.org/
 * 
 * This implementation is based on the code from BNichols available at:
 * http://processing.org/discourse/yabb2/YaBB.pl?num=1142535442
 */

package BeamCalc;

import java.io.File;

import processing.core.PApplet;
import processing.core.PVector;

/**
 * This class imports parses DXF files and converts the line
 * objects into element objects understandable by the solver
 * @author Enrique Ramos
 *
 */
public class DXFImporter {

	PApplet p;
	Structure structure;
	
	String[] dxf;
	String[] vport;
	String[] entity;
	String[][] ent;
	String[][] code;

	float viewcenterX = 0;
	float viewcenterY = 0;
	float viewWidthX = 1;
	float viewWidthY = 1;
	float aspectRatio = 1;
	float zoom = 1;

	DXFImporter(PApplet p, Structure struc) {
		this.structure = struc;
		this.p = p;
	}

	void DXFImport(File file) {
		readDXF(file);

		Float f = new Float(vport[0]);
		viewcenterX = f.floatValue();
		f = new Float(vport[2]);
		viewcenterY = f.floatValue();
		f = new Float(vport[28]);
		viewWidthY = f.floatValue();
		f = new Float(vport[30]);
		aspectRatio = f.floatValue();
		viewWidthX = viewWidthY * aspectRatio;
		// zoom = fixedheight / viewWidthY;
		
		importDxf();
		structure.draw(p);
	}

	void importDxf() {
		
		structure.reset();

	  for (int i = 1; i < ent.length; i++) {
	    if (ent[i][1].contains("LINE")) {

	    	Float f1 = new Float(ent[i][13]);
	    	Float f2 = new Float(ent[i][15]);
	    	Float f3 = new Float(ent[i][17]);
	    	Float f4 = new Float(ent[i][19]);
	    	Float f5 = new Float(ent[i][21]);
	    	Float f6 = new Float(ent[i][23]);
	    	
	    	PVector node1 = new PVector(f1.floatValue(), f2.floatValue(), f3.floatValue());
	    	PVector node2 = new PVector(f4.floatValue(), f5.floatValue(), f6.floatValue());
	    	
	    	
	    	
	    	if(Math.abs(f3.floatValue())<100){
	    		structure.addBeam(node1, node2, 3, 0, 2);
	    	}else if(Math.abs(f6.floatValue())<100){
	    		structure.addBeam(node1, node2, 0, 3, 2);
	    	}else{
	    		structure.addBeam(node1, node2, 0, 0, 2);
	    	}
	    	
	    }

	  }
	}

	void readDXF(File f) {
		//File f = new File(p.dataPath("envelope01.dxf"));
		dxf = PApplet.loadStrings(f);
		vport = cutSection(dxf, "VPORT", "ENDTAB");
		vport = cutSection(vport, " 12", " 43");
		dxf = cutSection(dxf, "ENTITIES", "ENDSEC");

		int numEntities = 0;
		for (int i = 0; i < dxf.length; i++) {
			if (dxf[i].contains("  0")) {
				dxf[i] = "ENTITY";
				numEntities++;
			}
		}
		String joindxf;
		joindxf = PApplet.join(dxf, "~");

		entity = PApplet.split(joindxf, "ENTITY");
		ent = new String[numEntities + 1][];
		for (int i = 0; i <= numEntities; i++) {
			ent[i] = PApplet.split(entity[i], "~");
		}
	}

	String[] cutSection(String[] dxfs, String startcut, String endcut) {
		int cutS = -1;
		for (int i = 0; i < dxfs.length; i++) {
			if (dxfs[i].contains(startcut)) {
				cutS = i;
			}
		}
		if (cutS == -1) {
			PApplet.println("SECTION " + startcut + " NOT FOUND.");
		}
		dxfs = PApplet.subset(dxfs, cutS + 1);

		int cutF = -1;
		for (int i = 0; i < dxfs.length; i++) {
			if (dxfs[i].contains(endcut)) {
				cutF = i;
				break;
			}
		}
		if (cutF == -1) {
			PApplet.println("SECTION NOT TERMINATED at " + endcut + ".");
		}
		return PApplet.subset(dxfs, 0, cutF - 1);
	}
}
