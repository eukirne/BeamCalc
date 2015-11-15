/**
 * (./) OpenFileAction.java v0.1 05/09/2011
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
import java.io.File;

import javax.swing.AbstractAction;
import javax.swing.JFileChooser;
import javax.swing.JFrame;

/**
 * This action creates and shows a modal open-file dialog.
 * @author Enrique Ramos
 *
 */
public class OpenFileAction extends AbstractAction {
    JFrame frame;
    JFileChooser chooser;

    OpenFileAction(JFrame frame, JFileChooser chooser) {
        super("Open...");
        this.chooser = chooser;
        this.frame = frame;
    }

    public void actionPerformed(ActionEvent evt) {
        // Show dialog; this method does not return until dialog is closed
        chooser.showOpenDialog(frame);

        // Get the selected file
        File file = chooser.getSelectedFile();
    }
};

