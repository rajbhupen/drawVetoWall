#include <iostream>
#include <vector>
#include <cmath>
#include "TCanvas.h"
#include "TBox.h"

//Just top veto wall picture
void drawMuonVetoWallBasic_v0() {
  
  double scintillatorWidth = 5.0;       // cm
  double scintillatorLength = 450.0;   // cm
  double scintillatorHeight = 1.0;     // cm
  double gapBetweenScintillators = 0.1; // cm
  double gapBetweenTiles = 0.2;        // cm
  double layerSpacing = 1.0;           // cm
  double fiberRadius = 2*0.07;  // Radius of fiber holes in cm
  double holeOffsetX = 1.25;   // Horizontal offset from EPS center
  double holeOffsetY = 0.5;   // Vertical offset from top and bottom
  int nEPS = 8;
  int nLayers = 4;
  int nTiles = 2;
  double layerStaggering = 1.25; // cm
  double layerHeights[4] = {2.0, 1.0, 2.0, 1.0};

  double totalWidth = 2 * (nEPS * scintillatorWidth + (nEPS - 1) * gapBetweenScintillators) + gapBetweenTiles + 3 * layerStaggering;
  double totalHeight = 2 * (1.0 + 1.0) + 2 * (1.0 + 2.0);

  // Arrays to store scintillator positions
  double xLeft[100];
  double xRight[100];
  double yBottom[100];
  double yTop[100];

  int nMuons = 5;
 
  int colorcode[10] = {1,2,3,4,6,7,8,9,46,47};


  TCanvas *c1 = new TCanvas("c1", "Muon Veto Wall Geometry", 30*totalWidth, 30*totalHeight);
  c1->Range(-5.0, -2.0, totalWidth + 5.0, totalHeight + 2.0);



  // Start y-position for the first layer
  double startY = 0.0;
  int scintIndex = 0;
  for (int layer = 0; layer < nLayers; ++layer) {
   
    double currentHeight = layerHeights[layer];
    double startX = (layer + 1) * layerStaggering;

    for (int tile = 0; tile < nTiles; ++tile) {
      for (int scint = 0; scint < nEPS; ++scint) {
	xLeft[scintIndex] = startX;
	xRight[scintIndex] = startX + scintillatorWidth;
	yBottom[scintIndex] = startY;
	yTop[scintIndex] = startY + currentHeight;

	TBox *scintillator = new TBox(xLeft[scintIndex], yBottom[scintIndex],
				      xRight[scintIndex], yTop[scintIndex]);
	scintillator->SetFillColor(kGreen - 8); 
	scintillator->Draw();

	for (int i = 0; i < 2; ++i) {
	  double holeX = xLeft[scintIndex] + scintillatorWidth/2.;
	  if(i==0) holeX-=holeOffsetX;
	  else holeX+=holeOffsetX;
		  
	  double holeY = yBottom[scintIndex] + currentHeight / 2.0;  
	  TEllipse *fiberHole = new TEllipse(holeX, holeY, fiberRadius);
	  fiberHole->SetFillColor(kWhite);
	  fiberHole->SetLineColor(kBlack);
	  fiberHole->Draw();
	}
	startX += scintillatorWidth + gapBetweenScintillators;
	scintIndex++;
      }
      startX += gapBetweenTiles;
    }
    startY += currentHeight + layerSpacing;
    
  }

    
  // Save and display the canvas
  c1->Update();
  c1->SaveAs("MuonVetoWall.png");
}

/*

  void drawMuonVetoWallBasic_v1() {
  double scintillatorWidth = 5.0;       // cm
  double scintillatorLength = 450.0;   // cm
  double scintillatorHeight = 1.0;     // cm
  double gapBetweenScintillators = 0.1; // cm
  double gapBetweenTiles = 0.2;        // cm
  double layerSpacing = 1.0;           // cm
  double fiberRadius = 2*0.07;  // Radius of fiber holes in cm
  double holeOffsetX = 1.25;   // Horizontal offset from EPS center
  double holeOffsetY = 0.5;   // Vertical offset from top and bottom
  int nEPS = 8;
  int nLayers = 4;
  int nTiles = 2;
  double layerStaggering = 1.25; // cm
  double layerHeights[4] = {2.0, 1.0, 2.0, 1.0};

  double totalWidth = 2 * (nEPS * scintillatorWidth + (nEPS - 1) * gapBetweenScintillators) + gapBetweenTiles + 3 * layerStaggering;
  double totalHeight = 2 * (1.0 + 1.0) + 2 * (1.0 + 2.0);

  // Arrays to store scintillator positions
  double xLeft[100];
  double xRight[100];
  double yBottom[100];
  double yTop[100];

  int nMuons = 5;
 
  int colorcode[10] = {1,2,3,4,6,7,8,9,46,47};


  TCanvas *c1 = new TCanvas("c1", "Muon Veto Wall Geometry", 30*totalWidth, 30*totalHeight);
  c1->Range(-5.0, -2.0, totalWidth + 5.0, totalHeight + 2.0);



  // Start y-position for the first layer
  double startY = 0.0;

  // Fill scintillator positions and draw wall
  int scintIndex = 0;
  for (int layer = 0; layer < nLayers; ++layer) {
   
  double currentHeight = layerHeights[layer];
  double startX = (layer + 1) * layerStaggering;

  for (int tile = 0; tile < nTiles; ++tile) {
  for (int scint = 0; scint < nEPS; ++scint) {
  xLeft[scintIndex] = startX;
  xRight[scintIndex] = startX + scintillatorWidth;
  yBottom[scintIndex] = startY;
  yTop[scintIndex] = startY + currentHeight;

  TBox *scintillator = new TBox(xLeft[scintIndex], yBottom[scintIndex],
  xRight[scintIndex], yTop[scintIndex]);
  scintillator->SetFillColor(kGreen - 8); // Default color
  scintillator->Draw();
  for (int i = 0; i < 2; ++i) {
  double holeX = xLeft[scintIndex] + scintillatorWidth/2.;
  if(i==0) holeX-=holeOffsetX;
  else holeX+=holeOffsetX;
		  
  double holeY = yBottom[scintIndex] + currentHeight / 2.0;  // Center vertically
  TEllipse *fiberHole = new TEllipse(holeX, holeY, fiberRadius);
  fiberHole->SetFillColor(kWhite);
  fiberHole->SetLineColor(kBlack);
  fiberHole->Draw();
  }
  startX += scintillatorWidth + gapBetweenScintillators;
  scintIndex++;
  }
  startX += gapBetweenTiles;
  }
  startY += currentHeight + layerSpacing;
    
  }


  // Generate and draw 5 random muon trajectories
  for (int i = 0; i < nMuons; ++i) {
  // Random starting position (x, y) in the bottom layer
  double xStart = gRandom->Uniform(0.0,totalWidth);
  
      
  double yStart = -1.0;                 // Fixed y position at bottom

  // Random direction (dx, dy)
  double dx = 2*gRandom->Uniform()-1; // Random dx from -1 to 1
  double dy = 1.0;                               // Fixed upward direction

  // Calculate intersections with scintillators
     
  TLine *muonLine = new TLine(xStart, yStart, xStart + 20 * dx, yStart + 20 * dy);
  muonLine->SetLineColor(colorcode[i]);
  muonLine->SetLineWidth(2);
      

  for (int j = 0; j < scintIndex; ++j) {
  cout<<j<<endl;
  bool isHit = false;

  if (dx != 0) {
  double tLeft = (xLeft[j] - xStart) / dx;
  double yLeft = yStart + tLeft * dy;
  if (tLeft >= 0 && yLeft >= yBottom[j] && yLeft <= yTop[j]) isHit = true;

  double tRight = (xRight[j] - xStart) / dx;
  double yRight = yStart + tRight * dy;
  if (tRight >= 0 && yRight >= yBottom[j] && yRight <= yTop[j]) isHit = true;
  }

  if (dy != 0) {
  double tBottom = (yBottom[j] - yStart) / dy;
  double xBottom = xStart + tBottom * dx;
  if (tBottom >= 0 && xBottom >= xLeft[j] && xBottom <= xRight[j]) isHit = true;

  double tTop = (yTop[j] - yStart) / dy;
  double xTop = xStart + tTop * dx;
  if (tTop >= 0 && xTop >= xLeft[j] && xTop <= xRight[j]) isHit = true;
  }

  // Color hit scintillator
  if (isHit) {
  muonLine->Draw();
  TBox *hitScint = new TBox(xLeft[j], yBottom[j], xRight[j], yTop[j]);
  hitScint->SetFillColor(colorcode[i]);
  //hitScint->SetFillStyle(3004); // Transparent fill
  hitScint->Draw();
  }
  }
  }
    
  // Save and display the canvas
  c1->Update();
  c1->SaveAs("MuonVetoWallWithPatterns.png");
  }
*/
void drawMuonVetoWallBasic_v2() {

  double scintillatorWidth = 5.0;       // cm
  double scintillatorLength = 450.0;   // cm
  double scintillatorHeight = 1.0;     // cm
  double gapBetweenScintillators = 0.1; // cm
  double gapBetweenTiles = 0.2;        // cm
  double layerSpacing = 1.0;           // cm
  double fiberRadius = 2*0.07;    // Radius of fiber holes in cm
  double holeOffsetX = 1.25;           // Horizontal offset from EPS center
  double holeOffsetY = 0.5;            // Vertical offset from top and bottom
  int nEPS = 8;
  int nLayers = 4;
  int nTiles = 2;
  double layerStaggering = 1.25; // cm
  double layerHeights[4] = {2.0, 1.0, 2.0, 1.0};

  double totalWidth = 2 * (nEPS * scintillatorWidth + (nEPS - 1) * gapBetweenScintillators) + gapBetweenTiles + 3 * layerStaggering;
  double totalHeight = 2 * (1.0 + 1.0) + 2 * (1.0 + 2.0);

  // Arrays to store scintillator positions
  double xLeft[100];
  double xRight[100];
  double yBottom[100];
  double yTop[100];

  int nMuons = 5;
  int colorcode[10] = {1, 2, 3, 4, 6, 7, 8, 9, 46, 47};

  TCanvas *c1 = new TCanvas("c1", "Muon Veto Wall Geometry", 30 * totalWidth, 30 * totalHeight);
  c1->Range(-5.0, -2.0, totalWidth + 5.0, totalHeight + 2.0);

  //Draw X and Y axis
  // TLine* xaxis = new TLine(0.0,0.0,totalWidth,0.0);
  // xaxis->Draw();
  // TLine* yaxis = new TLine(0.0,0.0,0.0,totalHeight);
  // yaxis->Draw();

  // Start y-position for the first layer
  double startY = 0.0;

  // Fill scintillator positions and draw wall
  int scintIndex = 0;
  for (int layer = 0; layer < nLayers; ++layer) {
    double currentHeight = layerHeights[layer];
    double startX = (layer) * layerStaggering;

    for (int tile = 0; tile < nTiles; ++tile) {
      for (int scint = 0; scint < nEPS; ++scint) {
	xLeft[scintIndex] = startX;
	xRight[scintIndex] = startX + scintillatorWidth;
	yBottom[scintIndex] = startY;
	yTop[scintIndex] = startY + currentHeight;

	TBox *scintillator = new TBox(xLeft[scintIndex], yBottom[scintIndex],
				      xRight[scintIndex], yTop[scintIndex]);
	scintillator->SetFillColor(kGreen - 8); // Default color
	scintillator->Draw();
	for (int i = 0; i < 2; ++i) {
	  double holeX = xLeft[scintIndex] + scintillatorWidth/2.;
	  if(i==0) holeX-=holeOffsetX;
	  else holeX+=holeOffsetX;
		  
	  double holeY = yBottom[scintIndex] + currentHeight / 2.0;  // Center vertically
	  TEllipse *fiberHole = new TEllipse(holeX, holeY, fiberRadius);
	  fiberHole->SetFillColor(kWhite);
	  fiberHole->SetLineColor(kBlack);
	  fiberHole->Draw();
	}
	startX += scintillatorWidth + gapBetweenScintillators;
	scintIndex++;
      }
      startX += gapBetweenTiles;
    }
    startY += currentHeight + layerSpacing;
  }

    
  // Generate and draw 1 random muon trajectory
  for (int i = 0; i < 5; ++i) {
    // Random starting position (x, y) in the bottom layer
    double xStart = gRandom->Uniform(0.0, totalWidth);
    double yStart = -1.0; // Fixed y position at bottom
    std::cout << "Muon starting position: xStart = " << xStart << ", yStart = " << yStart << std::endl;

    // Random direction (dx, dy)
    double dx = 2 * gRandom->Uniform() - 1; // Random dx from -1 to 1
    double dy = 1.0;                        // Fixed upward direction

    if(i==0){
      xStart = 5.05 + layerStaggering;//
      yStart = -1.0;
      dx = 0;
      dy = 1.0;      
    }

    double norm = sqrt(dx * dx + dy * dy);
    dx /= norm;
    dy /= norm;
    std::cout << "Muon direction: dx = " << dx << ", dy = " << dy << " (normalized)" << std::endl;

    
    // Calculate intersections with scintillators
    TLine *muonLine = new TLine(xStart, yStart, xStart + 20 * dx, yStart + 20 * dy);
    muonLine->SetLineColor(colorcode[i]);
    muonLine->SetLineWidth(2);

    for (int j = 0; j < scintIndex; ++j) {
      bool isHit = false;
      std::cout << "Checking scintillator " << j << " edges:" << std::endl;

      // Define edges of the scintillator
      double edges[4][4] = {
	{xLeft[j], yBottom[j], xRight[j], yBottom[j]}, // Bottom edge
	{xRight[j], yBottom[j], xRight[j], yTop[j]},   // Right edge
	{xRight[j], yTop[j], xLeft[j], yTop[j]},       // Top edge
	{xLeft[j], yTop[j], xLeft[j], yBottom[j]}      // Left edge
      };

      // Check each edge
      for (int edge = 0; edge < 4; ++edge) {
	double x1 = edges[edge][0], y1 = edges[edge][1];
	double x2 = edges[edge][2], y2 = edges[edge][3];

	// Edge vector
	double ex = x2 - x1;
	double ey = y2 - y1;

	// Edge length
	double edgeLength = sqrt(ex * ex + ey * ey);

	// Normalize edge vector
	ex /= edgeLength;
	ey /= edgeLength;

	std::cout << "  Edge " << edge << ": (" << x1 << ", " << y1 << ") to (" << x2 << ", " << y2 << ")" << std::endl;

	// Solve for t and s
	double det = dx * ey - dy * ex;
	if (fabs(det) < 1e-9) {
	  std::cout << "    Parallel line detected, skipping edge." << std::endl;
	  continue; // Parallel, no intersection
	}

	double t = ((x1 - xStart) * ey - (y1 - yStart) * ex) / det;
	double s = ((xStart - x1) * dy - (yStart - y1) * dx) / det;

	std::cout << "    Intersection parameters: t = " << t << ", s = " << s << std::endl;

	// Estimate the intersection point (xi, yi)
	if (t >= 0) {
	  double xi = xStart + t * dx;
	  double yi = yStart + t * dy;

	  // Check if the intersection point (xi, yi) lies within the boundaries of the scintillator
	  bool withinBounds = (xi >= xLeft[j] && xi <= xRight[j] && yi >= yBottom[j] && yi <= yTop[j]);

	  if (withinBounds) {
	    std::cout << "    Intersection point: (" << xi << ", " << yi << ") is within the scintillator boundary." << std::endl;
	    isHit = true;
	    break;
	  } else {
	    std::cout << "    Intersection point: (" << xi << ", " << yi << ") is outside the scintillator boundary." << std::endl;
	  }
	}
      }

      // Visualize the hit
      if (isHit) {
	std::cout << "  Muon hits scintillator " << j << "." << std::endl;
	muonLine->Draw();
	TBox *hitScint = new TBox(xLeft[j], yBottom[j], xRight[j], yTop[j]);
	hitScint->SetFillColor(colorcode[i]); // Set hit color
	hitScint->Draw();
      } else {
	std::cout << "  Muon does not hit scintillator " << j << "." << std::endl;
      }
    }
  }
 
  // Save and display the canvas
  c1->Update();
  c1->SaveAs("MuonVetoWallWithPatterns.png");
}
