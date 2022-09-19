#include <iostream>       // std::cout
#include <string>         // std::string

#include <TMath.h>
#include <TRandom3.h>

#include <TH1.h>
#include <TF1.h>

#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>


void testPEwithNoise(int nEvents=100000, float noiseLevel=1, float avPEmult=1, float PEampl=5, float relPEwidth=0.25)
{ /* parameters: 
	 number of generated events, 
	 noise around pedestal (gaussian distribution, a.u.), 
	 average single PE pulse multiplicity (poisson distibution), 
	 (average) amplitude of injected pulse (a.u.), 
	 relative width of pulse amplitude distribution (gaussian distribution)
    */
  
  // each event consists of N time samples
  const int nTimeSamp = 100;
  int nSampInPed=nTimeSamp/2;   // only use this many samples for per-event pedestal calculation

  // --- modelling ---
  // background: each channel has a pedestal, and some noise around it --> sample a gaussian(pedestal, noiseLevel), put in a vector
  float pedestal=40;
  TF1 *fBkg=new TF1("fBkg","gaus(0)", -200, 200);
  fBkg->SetParameters(1000, pedestal, noiseLevel);
  // injected pulse: sample poisson to get a random pulse multiplicity per event around some average 
  TF1 *fPE=new TF1("fPE", "TMath::Poisson(x, [0])", -100, 100);
  fPE->SetParameter(0, avPEmult);
  // sample a gaussian around the average injected pulse amplitude. width defaults to 25%.
  // note that the pulse is always the same shape: two consecutive time samples, with the second at 40% amplitude compared to the first 
  TF1 *fPulse=new TF1("fPulse", "gaus", 0, 100);
  fPulse->SetParameters(1, PEampl, relPEwidth*PEampl);


  // --- setup histograms ---
  int binsPerInt=2; //binning resolution
  // amplitude sums
  // the relevant range is known from the pedestal and expected noise (which on the high side contains extra pulses) 
  int sumMinBin=(pedestal-4*noiseLevel)*nTimeSamp;
  int sumMaxBin=(pedestal+TMath::Max((int)(avPEmult*PEampl), 4)*noiseLevel)*nTimeSamp;
  int sumNBins = (sumMaxBin-sumMinBin)*binsPerInt;
  TH1F * hSum = new TH1F("hSum", "Simple sum histogram;#Sigma(amplitudes);Events", sumNBins, sumMinBin, sumMaxBin);
  // pedestal subtracted amplitude sums
  int subSumMinBin=(-TMath::Max((int)(avPEmult*PEampl), 4)*noiseLevel)*nTimeSamp;
  int subSumMaxBin=( TMath::Max((int)(avPEmult*PEampl), 4)*noiseLevel)*nTimeSamp;
  int subSumNBins = (subSumMaxBin-subSumMinBin)*binsPerInt;
  TH1F * hSubtractedSum = new TH1F("hSubtractedSum", "Pedestal subtracted sum histogram;#Sigma(pedestal subtracted amplitudes);Events", subSumNBins, subSumMinBin, subSumMaxBin); 
  TH1F * hMedSubtractedSum = new TH1F("hMedSubtractedSum", "Median range average pedestal subtracted sum histogram;#Sigma(pedestal subtracted amplitudes);Events", subSumNBins, subSumMinBin, subSumMaxBin); 
  // "truth"; the actual added charge in the event
  TH1F * hAddedCharge = new TH1F("hAddedCharge", "Added total pulse amplitude histogram;#Sigma(added amplitudes);Events", subSumNBins, subSumMinBin, subSumMaxBin);
  // cross check plot: time sample chosen for pulse injection; should tend to uniform 
  TH1F * hTimeSample = new TH1F("hTimeSample", "Start time sample chosen for PE pulses", nTimeSamp+1, 0, nTimeSamp+1);
  // allow 10 event displays for visual inspection
  TH1F * hEvent[10];
  float addedCharge[10]={0}; // for book keeping and printout

  // --- plotting ---
  // for printouts on plots
  TLatex * latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.04);
  
  TCanvas * c1 = new TCanvas("c1", "Single events canvas", 1200, 1200);
  c1->Divide(5,2);

  // used for pulling the time sample of injected pulse(s)
  TRandom3 * r = new TRandom3();

  //loop over events 
  for (int iE = 0; iE<nEvents; iE++) {
	vector <float> samples; // amplitude in each time sample
	float add=0;            // keep track of added pulse amplitudes total
	float sum = 0;          // sum of charge in the event
	float ped=0;            // calculated pedestal in the event

	if (iE < 10) //initialize event display for 10 first events
	  hEvent[iE] = new TH1F(Form("h%i", iE), Form("Readout in event %i;Time sample;Readout, ev %i", iE, iE), nTimeSamp+1,0, nTimeSamp+1);

	//loop over time samples 
	for (int iT = 0; iT<nTimeSamp; iT++) {
	  float val = fBkg->GetRandom(); //baseline value in the time sample
	  samples.push_back(val);
	  sum += val;
	  if (iT < nSampInPed) // only use some initial portion for average pedestal 
		ped += val;
	  if (iE < 10) 
		hEvent[iE]->Fill(iT+1, val);
	}//over time samples, setting up baseline
	
	// sample a poisson to get the number of PEs to put in the event, on top of the pedestal level readout.
	int nPEs=fPE->GetRandom();
	// that many times, pick a time sample t from a uniform distribution, smear ampl by gaussian, put 100% in t and 40% in t+1
	for (int iP=0; iP<nPEs; iP++) {
	  int ts = r->Uniform(0, nTimeSamp-2);
	  float amp=fPulse->GetRandom();
	  // add that to the corresponding vector index
	  samples[ts]=samples[ts]+1.4*amp;
	  // already truncated sampled range before the last time sample so we are sure to not go outside array; less than natural but sufficient
	  if (iE < 10)  {
		hEvent[iE]->Fill(ts+1, amp);
		hEvent[iE]->Fill(ts+1+1, 0.4*amp);
		addedCharge[iE]+=1.4*amp;
	  }
	  if (ts+1 < nSampInPed)
		ped+=0.4*amp;
	  sum+=1.4*amp;
	  add+=1.4*amp;
	  if (ts < nSampInPed)
		ped+=1.4*amp; 

	  hTimeSample->Fill(ts);
	}//over number of added pulses

	hAddedCharge->Fill(add);
	// take the (early) average --> this is the pedestal.
	ped/=nSampInPed;
	//alternate pedestal definition: average over the median two quantiles (exclude lowest and highest to reduce signal bias)
	std::sort(samples.begin(), samples.end());
    float medianPed=0;
	// use the same number of samples for the average-over-median method as for the average pedestal method  
	int pedOffset=nTimeSamp/2-nSampInPed/2; // center the interval 
	for (int i = pedOffset; i < pedOffset+nSampInPed ; i++) {
      medianPed+=samples[i];
	}
    medianPed/=nSampInPed;

	/* 
	   here it is probably useful to histogram the pedestals from the different methods 
	   (and maybe from a classical simple median method) to check performance, 
	   also histogram the injected pulse amplitude 
	*/
	
	// get the sum(vector) with and without subtracting the pedestal
	// histogram both.
	hSum->Fill(sum);	
	hSubtractedSum->Fill(sum-ped*nTimeSamp);	
	hMedSubtractedSum->Fill(sum-medianPed*nTimeSamp);	
	
	if (iE < 10) {
	  c1->cd(iE+1);
	  hEvent[iE]->Draw("hist");
	  latex->DrawLatex(0.2, 0.3,  Form("Added amplitude: %.2f", addedCharge[iE]) ); 
	  latex->DrawLatex(0.2, 0.25, Form("Found average pedestal : %.2f", ped) ); 
	  latex->DrawLatex(0.2, 0.2,  Form("Found median pedestal : %.2f", medianPed) ); 
	}
	
  }//over events

  
  // draw
  TCanvas * c4 = new TCanvas("c4", "PE pulse time sample canvas", 600, 500);
  hTimeSample->Draw("hist");

  TCanvas * c5 = new TCanvas("c5", "Subtracted sum canvas", 600, 500);
  c5->SetLogy();
  hSubtractedSum->Draw("hist");
  hMedSubtractedSum->SetLineColor(kAzure-9);
  hMedSubtractedSum->Draw("hist same");
  // "truth"
  hAddedCharge->SetLineColor(kGreen+2);
  hAddedCharge->SetLineStyle(2);
  hAddedCharge->Draw("hist same");
  cout << "\nFitting pedestal-subtracted sum of amplitudes" << endl;
  TF1 * fG2 = new TF1("fG2", "gaus", hSubtractedSum->GetXaxis()->GetXmin(),  hSubtractedSum->GetXaxis()->GetXmax());
  hSubtractedSum->Fit(fG2, "NR", "", -nTimeSamp*2*noiseLevel,nTimeSamp*2*noiseLevel);
  fG2->SetLineColor(kRed+1);
  fG2->SetLineWidth(3);
  fG2->Draw("same");


  TLegend * leg = new TLegend(0.7, 0.67, 0.95, 0.95);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  leg->AddEntry(hSubtractedSum, "Av. ped", "L");
  leg->AddEntry(fG2, "Gauss fit to av. ped.", "L");
  leg->AddEntry(hMedSubtractedSum, "Av. median ped", "L");
  leg->AddEntry(hAddedCharge, "Added pulse sum", "L");
  
  leg->Draw();

  
  TCanvas * c3 = new TCanvas("c3", "Subtracted sum canvas", 600, 500);
  c3->SetLogy();
  hSubtractedSum->Draw("hist");
  fG2->Draw("same");

   TCanvas * c2 = new TCanvas("c2", "Sum canvas", 600, 500);
  c2->SetLogy();
  hSum->Draw("hist");
  cout << "\nFitting sum of amplitudes" << endl;

  TF1 * fG1 = new TF1("fG1", "gaus", hSum->GetXaxis()->GetXmin(),  hSum->GetXaxis()->GetXmax());
hSum->Fit(fG1, "NR", "",hSum->GetXaxis()->GetXmin(), nTimeSamp*(pedestal+2*noiseLevel));
  fG1->SetLineColor(kCyan);
  fG1->SetLineWidth(5);
  fG1->Draw("same");


// done.


  
}

  
