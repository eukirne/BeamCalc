/**
 * (./) Toolbox.java v0.1 05/09/2011
 * @author Enrique Ramos Melgar
 * http://www.esc-studio.com
 *
 * THIS LIBRARY IS RELEASED UNDER A CREATIVE COMMONS ATTRIBUTION 3.0 LICENSE
 * http://creativecommons.org/licenses/by/3.0/
 * http://www.processing.org/
 * 
 */

package BeamCalc;


import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.*;
import javax.swing.event.*;
/**
 * This class creates a java JFrame containing a Processing application
 * containing all the controls of the solver.
 * @author Enrique Ramos
 */
public class Toolbox {
	Structure structure;
	JSlider slider, slider2, slider3, slider4, slider5, slider6, slider7;
	JLabel label, label2, label3, label4, label5, label6, label7;
	JButton button1, button2, button3, button4, button5, button6, button7, button8,button9, button10;
	JFrame frame;
	
	JFileChooser chooser;
	
	int xPos = 10;
	int yPos = 10;
	int buttonW = 90;

	public Toolbox(final Structure structure) {
		this.structure = structure;
		frame = new JFrame("Parameters");
		
		slider = new JSlider();
		slider.setMinimum(0);
		slider.setMaximum(500);
		slider.setValue(20);
		slider.addChangeListener(new MyChangeAction());
		
		slider2 = new JSlider();
		slider2.setMinimum(1);
		slider2.setMaximum(100);
		slider2.setValue(10);
		slider2.addChangeListener(new MyChangeAction());
		slider2.setBounds(100, 100, 100, 10);
		
		// Step
		slider3 = new JSlider();
		slider3.setMinimum(0);
		slider3.setMaximum(1000);
		slider3.setValue(400);
		slider3.addChangeListener(new MyChangeAction());
		slider3.setBounds(100, 100, 100, 10);
		
		slider4 = new JSlider();
		slider4.setMinimum(0);
		slider4.setMaximum(100);
		slider4.setValue(50);
		slider4.addChangeListener(new MyChangeAction());
		slider4.setBounds(100, 100, 100, 10);
		
		slider5 = new JSlider();
		slider5.setMinimum(10);
		slider5.setMaximum(500);
		slider5.setValue(80);
		slider5.addChangeListener(new MyChangeAction());
		slider5.setBounds(100, 100, 100, 10);
		
		slider6 = new JSlider();
		slider6.setMinimum(0);
		slider6.setMaximum(50);
		slider6.setValue(15);
		slider6.addChangeListener(new MyChangeAction());
		slider6.setBounds(100, 100, 100, 10);
		
		slider7 = new JSlider();
		slider7.setMinimum(0);
		slider7.setMaximum(100);
		slider7.setValue(0);
		slider7.addChangeListener(new MyChangeAction());
		slider7.setBounds(100, 100, 100, 10);
		
		label = new JLabel("Displacement Scale");
		label2 = new JLabel("Moment Scale");
		label3 = new JLabel("Step");
		label4 = new JLabel("Section Area");
		label5 = new JLabel("Number of Nodes");
		label6 = new JLabel("Flocking Distance");
		label7 = new JLabel("Wind Load");
		
		JPanel panel = new JPanel();
		panel.add(slider);
		panel.add(label);
		panel.add(slider2);
		panel.add(label2);
//		panel.add(slider7);
//		panel.add(label7);
//		panel.add(slider3);
//		panel.add(label3);
//		panel.add(slider4);
//		panel.add(label4);
//		panel.add(slider5);
//		panel.add(label5);
//		panel.add(slider6);
//		panel.add(label6);
		
		
		yPos = frame.getHeight()/2+120;
		
		button1 = new JButton("ImpDXF");
		button1.setSize(buttonW, 30);
		button1.setLocation(xPos,yPos);
		
		button6 = new JButton("Topology");
		button6.setSize(buttonW, 30);
		button6.setLocation(xPos+buttonW,yPos); yPos+=30;
		
		button8 = new JButton("Editing");
		button8.setSize(buttonW, 30);
		button8.setLocation(xPos,yPos);
		
		
		button7 = new JButton("Profile On");
		button7.setSize(buttonW, 30);
		button7.setLocation(xPos + buttonW,yPos);yPos+=30;
		
		button9 = new JButton("Graph");
		button9.setSize(buttonW, 30);
		button9.setLocation(xPos,yPos);
		
		button10 = new JButton("View All");
		button10.setSize(buttonW, 30);
		button10.setLocation(xPos+ buttonW,yPos); yPos+=30;
		
		button2 = new JButton("Disp");
		button2.setSize(buttonW, 30);
		button2.setLocation(xPos,yPos);

		button3 = new JButton("Moment");
		button3.setSize(buttonW, 30);
		button3.setLocation(xPos+ buttonW,yPos); yPos+=30;

		button4 = new JButton("Axial");
		button4.setSize(buttonW, 30);
		button4.setLocation(xPos,yPos);
		
		button5 = new JButton("Torsion");
		button5.setSize(buttonW, 30);
		button5.setLocation(xPos+ buttonW,yPos); yPos+=30;
		
		
		
		
		frame.add(button1);
		frame.add(button2);
		frame.add(button3);
		frame.add(button4);
		frame.add(button5);
		frame.add(button6);
		frame.add(button7);
		frame.add(button8);
		frame.add(button9);
		frame.add(button10);
//		
		frame.add(panel);
	    frame.setSize(200, structure.p.frame.getHeight()/2-5);
	    frame.setLocation(structure.p.frame.getLocation().x + structure.p.frame.getWidth() + 5,structure.p.frame.getLocation().y+structure.p.frame.getHeight()/2+5);
	    frame.setVisible(true);
	    
		// Add action listener to button
		button1.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				structure.analysing=false;
				structure.drawing = false;
				structure.iterations=1;
				JFrame fileFrame = new JFrame();
			    JFileChooser chooser = new JFileChooser();
				chooser.showOpenDialog(fileFrame);

		        // Get the selected file
		        File file = chooser.getSelectedFile();
		        if(file!=null){
		        	structure.nodes.removeAll(structure.nodes);
		        	structure.elms.removeAll(structure.elms);
		        	structure.dxf.DXFImport(file);
		        }
		        structure.analysing=true;
				structure.drawing = true;
				structure.p.frame.toFront();
			}
		});
		// Add action listener to button
		button2.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				structure.displayMode=4;
				structure.p.frame.toFront();
			}
		});
		// Add action listener to button
		button3.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				structure.displayMode=2;
				structure.p.frame.toFront();
			}
		});
		// Add action listener to button
		button4.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				structure.displayMode=1;
				structure.p.frame.toFront();
			}
		});
		// Add action listener to button
		button5.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				structure.displayMode=3;
				structure.p.frame.toFront();
			}
		});
		// Add action listener to button
		button6.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				structure.changeTopology();
				structure.p.frame.toFront();
			}
		});
		// Add action listener to button
		button7.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				if(structure.rendering){
					structure.render(false);
					button7.setText("Profile Off");
				}else{
					structure.render(true);
					button7.setText("Profile On");
				}
				structure.p.frame.toFront();
			}
		});
		// Add action listener to button
		button8.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				structure.displayMode=0;
				structure.p.frame.toFront();
			}
		});
		button9.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				structure.graph =!structure.graph;
				if(structure.graph){
					button9.setText("Graph On");
				}else{
					button9.setText("Graph Off");
				}
				structure.p.frame.toFront();
			}
		});
		button10.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				structure.focusAll();
				structure.p.frame.toFront();
			}
		});
	}
	
	
	
	float getDispScale(){
    	return slider.getValue();
    }
	void setGravity(float value){
		slider.setValue((int)value);
	}
	
	float getMomentScale(){
    	return (float)slider2.getValue()/10;
    }
	void setK(float value){
		slider2.setValue((int)value);
	}
	
	float getStep(){
    	return (float)slider3.getValue()/100000000;
    }
	void setStep(float value){
		slider3.setValue((int)(value));
	}
	
	float getArea(){
    	return (float)slider4.getValue()/100;
    }
	void setArea(float value){
		slider4.setValue((int)(value));
	}
	
	float getNumNodes(){
    	return (int)slider5.getValue();
    }
	void setNumNodes(float value){
		slider5.setValue((int)value);
	}
	
	float getFlockDist(){
    	return (float)slider6.getValue();
    }
	void setFlockDist(float value){
		slider6.setValue((int)value);
	}
	
	float getThreshold(){
    	return (float)slider7.getValue()/100;
    }
	void setThreshold(float value){
		slider7.setValue((int)(value*100));
	}

	public class MyChangeAction implements ChangeListener {
		public void stateChanged(ChangeEvent ce) {
			int value = slider.getValue();
			String str = Integer.toString(value);
			label.setText("Displacement Scale: " + str);
			value = slider2.getValue();
			str = Integer.toString(value);
			label2.setText("Moment Scale: " + str);
			value = slider3.getValue();
			str = Integer.toString(value);
			label3.setText("Step: " + "0.00" + str +" m");
			value = slider4.getValue();
			str = Integer.toString(value);
			label4.setText("Section Area: " + str);
			value = slider5.getValue();
			str = Integer.toString(value);
			label5.setText("Number of Nodes: " + str);
			value = slider6.getValue();
			str = Integer.toString(value);
			label6.setText("Flocking Distance: " + str);
			value = slider7.getValue();
			str = Integer.toString(value);
			label7.setText("Wind Load: " + str);
		}
	}
}