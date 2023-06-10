#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1F.h"


void getBeamTrajectory(float position=-880.)  {

  //known points 
  TGraphErrors * gPath = new TGraphErrors(3);
  TGraphErrors * gAngle = new TGraphErrors(3);

  //use the points along the trajectory used/assumed in simulation
  gPath->SetPoint(0, -880., -44.); 
  gPath->SetPoint(1, -700., -27.926);
  gPath->SetPoint(2, 0.,0.);

  //use the angles along the trajectory used/assumed in simulation
  gAngle->SetPoint(0, -880., 5.65);
  gAngle->SetPoint(1, -700., 4.5);
  gAngle->SetPoint(2, 0.,0.);

  //some plotting 
  gPath->Draw("ap");
  gPath->GetHistogram()->GetYaxis()->SetRangeUser(-50, 10);
  gPath->GetHistogram()->GetXaxis()->SetRangeUser(-1000, 10);
  gPath->DrawClone("p");
  gPath->SetTitle(";z [mm];beam trajectory x (mm, black) or #theta (#circ, red)");

  gAngle->SetMarkerColor(kRed);
  gAngle->SetLineColor(gAngle->GetMarkerColor());
  gAngle->SetMarkerStyle(kOpenSquare);
  gAngle->DrawClone("p same");

  //now fit to extract the beding radius 
  TF1 * fBeam = new TF1("fBeam", "-[0]+sqrt( [0]*[0] - (x*x) )", -1000., 100);
  fBeam->SetParameter(0, -8800); //some support needed; manually found -8787 and -8822 respectively from input points 0 and 1 
  gPath->Fit(fBeam);
  fBeam->DrawCopy("same");

  // empirically this looks like a straight line 
  TF1 * fAngle = new TF1("fAngle", "[0]+[1]*x", -1000., 100);
  gAngle->Fit(fAngle);
  fAngle->SetLineColor(gAngle->GetMarkerColor());
  fAngle->DrawCopy("same");

  //print the values at the input point used as argument
  std::cout << "From the fits, we should expect the best x at z = " << position << " to be " << fBeam->Eval(position) << " and the angle of incidence to be " << fAngle->Eval(position) << " degrees." << std::endl;

  //Done.
}
