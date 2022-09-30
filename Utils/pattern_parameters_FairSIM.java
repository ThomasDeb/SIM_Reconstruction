import java.util.List;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Scanner;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import org.fairsim.fiji.DisplayWrapper;
import org.fairsim.linalg.Vec2d;
import org.fairsim.linalg.Vec2d.Real;
import org.fairsim.sim_algorithm.OtfProvider;
import org.fairsim.sim_algorithm.SimAlgorithm;
import org.fairsim.sim_algorithm.SimParam;
import org.fairsim.sim_algorithm.SimUtils;
import org.fairsim.utils.ImageDisplay;
import org.fairsim.utils.Tool;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;


public class test{

	public static void main(String[] args) {


		// Path to image stack
		String path = "/Users/thomas/Documents/SIM_Reconstruction/Data_Arne/goodPixelSize/Nice/test_SIM026/test_SIM026_stack.tif";

		new ImageJ();
		ImagePlus imp = IJ.openImage(path);
		imp.show();

		// Hardware parameters
		int nrBands = 2;
		int nrAng = 3;
		int nrPha = 3;

		float background = 0;
		double pxlSize = 0.064;
		double emWavelen = 605;
		double otfNA = 1.49;	// NA
		double otfCorr = 0.1;	// compensation for the non-real OTF

		write_params(imp, nrBands, nrAng, nrPha, background, pxlSize, emWavelen, otfNA, otfCorr);
		imp.close();

	}

	private static void write_params(ImagePlus imp, int nrBands, int nrAng, int nrPha, float background, double pxlSize, double emWavelen, double otfNA, double otfCorr) {
		ImageStack iSt = imp.getStack();
		ImageDisplay.Factory idf = DisplayWrapper.getFactory();

		int imgSize = imp.getWidth()
		int nbTime = (int) (imp.getNSlices() / (nrAng * nrPha));
		int fitBand = 2;
		int visualFeedback = 0;
		int b = nrBands-1;

		double fitExclude = 0.6;		// freq. region (DC) to ignore when search for k0, in fraction of OTF cutoff

		boolean otfBeforeShift = false;	// if to apply the OTF before or after shifting the band


		// copy raw data, window, fft
		Vec2d.Cplx [][] rawImages = new Vec2d.Cplx[nrAng][nrPha];

		for (int t = 0; t<nbTime; t++) {
			for (int a = 0; a<nrAng; a++) {
				for (int p=0; p<nrPha; p++) {

					// get the current image as 16bit short
					short [] curImg = (short[])iSt.getProcessor( t*nrPha*nrAng + a*nrPha + p + 1  ).convertToShortProcessor().getPixels();

					// copy to a complex-value vector
					rawImages[a][p] = Vec2d.createCplx( imgSize, imgSize );
					rawImages[a][p].setFrom16bitPixels( curImg );

					// subtract background, window, fft
					double pc = SimUtils.subtractBackground( rawImages[a][p], background );

					SimUtils.fadeBorderCos( rawImages[a][p], 15);
					rawImages[a][p].fft2d(false);

					Tool.trace("fft'd input "+a+" "+p+", subtracted background, % pixels clipped: "+pc*100);
				}
			}
			// setup OTF and SIM parameters. Both of these could be loaded from an xml file!
			OtfProvider otf   = OtfProvider.fromEstimate( otfNA, emWavelen, otfCorr );
			SimParam simParam = SimParam.create( nrBands, nrAng, nrPha, imgSize, pxlSize, otf );

			// run the parameter fit
			SimAlgorithm.estimateParameters( simParam, rawImages, fitBand, fitExclude, idf, visualFeedback, null);

			for (int a = 0; a<nrAng; a++) {
				double [] k = simParam.dir(a).getPxPy(b);
				double sign = Math.signum(k[0]);
				simParam.dir(a).setPxPy(sign*k[0], sign*k[1]);
			}



			//	Now we write the parameters in a text file
			try {
				FileWriter fw = new FileWriter("/Users/thomas/Documents/SIM_Reconstruction/Data_Arne/goodPixelSize/Nice/test_SIM026/fairSIM.txt", true);
				BufferedWriter bw = new BufferedWriter(fw);
				PrintWriter pw = new PrintWriter(bw);
				// Write params in txt. Y is negative to account for fairSIM systematic coordinate change
				pw.println(simParam.dir(0).px(b)+" " + (simParam.dir(0).py(b)) +" "+ simParam.dir(0).getPhaOff());
				pw.println(simParam.dir(1).px(b)+" " + (simParam.dir(1).py(b))+" "+ simParam.dir(1).getPhaOff());
				pw.println(simParam.dir(2).px(b)+" " + (simParam.dir(2).py(b))+" "+ simParam.dir(2).getPhaOff());
				pw.flush();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				System.out.println("Data not appended");;
			} finally {System.out.println("Finally Statement"); }

			//	simParam.Dir.getPhOff;
		}
	}
}
