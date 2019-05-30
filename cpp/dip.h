#ifndef DIP_H
#define DIP_H

#include "var.h"

using namespace cv;


class DIP
{
  public:
	DIP(Mat receivedImage);
    static QImage mat2QImage(Mat& in);
	void processLaserStripeImage();
	void getImageSize();
	void getGrayImage();
	void getROI();
	void thinLaserStripeImage();
	void findFeaturePoints();
	double det(Point2i a, Point2i b);
	void generateInformation();
	void generateOutImage();
	void drawAsterisk(Mat& image, Point pt, Scalar& color, int thickness);



	

public:
	Mat image;
	int width;
	int height;
	Mat gray;
	Mat roi;
	static int ROIX;
 	static int ROIY;
	int roiHeight;
	int roiWidth;
	Mat filteredImage;
	Mat filteredImageDoubleType;
	Mat thinedImageInfo;
	Mat thinedImage;
	Point2i featurePointA;
	Point2i featurePointB;
	Point2i featurePointC;
	Point2i featurePointD;
	Point2i startPoint;
	Point2i endPoint;
	std::vector<Point2i> featurePoints;
	int offset;
	int weldSeamWidth;
	int weldSeamDepth;
	Mat outputImage;
	QImage out;
	std::vector<double>maxLoc;
	double slope;
	double maxSlope = 0;
	double minSlope = 0;
	int maxSlopeLoc = 0;
	int minSlopeLoc = 0;
	double t1 = 0, t2 = 0, t3 = 0, t4 = 0, b0 = 0, b1 = 0, yi = 0;
	std::vector<Point2i>centerLine;
	Point2i point;

};




#endif // DIP_H
