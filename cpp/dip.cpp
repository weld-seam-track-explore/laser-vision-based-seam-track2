/*****************************************************************************
Copyright: Guangdong Provincial Welding Technology Research Center
Author: Li Yangjin
Date: 2019-5-29
Description:DIP类对传感图像进行图像处理和特征点识别。其输入时Mat类型图形
，输出焊缝中心和焊缝偏差、焊缝深度和焊缝宽度。由于openCV和QT包含了不同的图像类型，DIP类
还提供了Mat和QImage格式转换的相关接口。

特征点提取算法:
	1.裁剪ROI
	2.滤波（包括图像滤波、弧光信息去除等）
	3.直线拟合（要包括细化）
	4.图像变换（前提是结构光条纹线大体保持水平）
	5.特征点提取

*****************************************************************************/
#include "dip.h"


int DIP::ROIX = 130;
int DIP::ROIY = 120;



/*****************************************************************************
Function:void DIP(Mat& recImg)
Description:对传入图形进行处理
Call:
Input:recImg-待处理的图像
Output:
Return:
Others:
*****************************************************************************/
DIP::DIP(Mat receivedImage)
{
	receivedImage.copyTo(image);
	processLaserStripeImage();

}



/*****************************************************************************
Function:DIP::processLaserStripeImage()
Description:根据焊缝的不同类型，采用不同的处理办法

	包括1. 裁剪ROI、2.滤波、3.直线拟合、4.图像变换、5.特征点提取

Input:
Output:
Return:
Others:
*****************************************************************************/


void DIP::processLaserStripeImage()
{
	getImageSize();
	getGrayImage();
	getROI();
	medianBlur(roi, filteredImage, 3);
	thinLaserStripeImage();
	findFeaturePoints();
	generateInformation();
	generateOutImage();

}






/*****************************************************************************
Function:void DIP::getImageSize()
Description:获得原图尺寸
Input:
Output:
Return:
Others:
*****************************************************************************/

void DIP::getImageSize()
{
	width = image.cols;
	height = image.rows;
}

/*****************************************************************************
Function:void DIP::getGrayImage()
Description:生成灰度图
Input:
Output:
Return:
Others:
*****************************************************************************/

void DIP::getGrayImage()
{
	gray = Mat::zeros(Size(image.cols, image.rows), CV_8UC1);
	cvtColor(image, gray, CV_BGR2GRAY);
}


/*****************************************************************************
Function:void DIP::getROI()
Description:裁剪原图,获得ROI区域，是为了降低运算量。根据投影自适应地裁剪原图，产生大小
	为250*130的ROI图。
Input:
Output:总体上，获得了roi,roiWidth，roiHeight
	这些内容。
Return:
Others:
*****************************************************************************/
void DIP::getROI()
{
	roiHeight = 130;
	roiWidth = 250;
	roi = gray(Range(ROIY, ROIY + roiHeight), Range(ROIX, ROIX + roiWidth));

}

/*****************************************************************************
Function:void DIP::thinLaserStripeImage()
Description:细化图像
Input:filteredImage
Output:thinnedImage，其为对filteredImage进行细化后的图像。
Return:
Others:
*****************************************************************************/

void DIP::thinLaserStripeImage()
{
	/*filteredImage.convertTo(filteredImageDoubleType, CV_64FC1);
	Mat assistA = Mat::ones(Size(roiHeight, 1), CV_64FC1);
	std::vector<double> assist;
	for (int i = 0; i < roiHeight; i++) assist.push_back(i);
	Mat assistB = Mat(assist);
	Mat assistBT = assistB.t();

	Mat num = assistBT * filteredImageDoubleType;
	Mat den = assistA * filteredImageDoubleType;
	thinedImageInfo = num / den;
	thinedImageInfo.convertTo(thinedImageInfo, CV_8UC1);
	thinedImage = Mat::zeros(Size(roiWidth, roiHeight), CV_8UC1);

	uchar* ptr = thinedImageInfo.ptr<uchar>(0);
	for (int i = 0; i < roiWidth; i++)
	{
		if (ptr[i])
			thinedImage.at<uchar>(ptr[i], i) = 255;
	}*/
	
	for (int i = 0; i < roiWidth; i++)
	{
		uint yMaxValue = 0;
		for (int j = 0; j < roiHeight; j++)
		{
			uchar* ptr = filteredImage.ptr<uchar>(j);
			if (yMaxValue < ptr[i]) yMaxValue = ptr[i];

		}

		int locSum = 0;
		int locNum = 0;
		double locMean = 0;
		for (int j = 0; j < roiHeight; j++)

		{
			uchar* ptr = filteredImage.ptr<uchar>(j);
			if (ptr[i] >= yMaxValue)
			{
				locSum += j;
				locNum++;
			}

		}
		locMean = locSum / locNum;
		maxLoc.push_back(locMean);


	}

	thinedImageInfo = (Mat)maxLoc;
	thinedImageInfo.convertTo(thinedImageInfo, CV_8UC1);
	thinedImage = Mat::zeros(Size(roiWidth, roiHeight), CV_8UC1);
	uchar* ptr = thinedImageInfo.ptr<uchar>(0);
	for (int i = 0; i < roiWidth; i++)
	{
		if (ptr[i])
		{
			thinedImage.at<uchar>(ptr[i], i) = 255;
		}
	}




}


/*****************************************************************************
Function:void DIP::findFeaturePoints()
Description:找焊缝特征点。利用斜率法算法从thinnedImageInfo获得焊缝特征点，根据
特征点将激光条纹分成3组，每组用最小二乘法拟合中心线
Input:
Output:
Return:返回找到的焊缝特征点，对于对接焊缝，返回两个特征点
Others:
*****************************************************************************/
void DIP::findFeaturePoints()
{
	//int width = thinedImageInfo.cols;
	//Point2i startPoint;
	//Point2i endPoint;
	//Point2i vertexPoint;

	//uchar* ptr = thinedImageInfo.ptr<uchar>();
	//uint yvalue = ptr[0];
	//for (int i = 0; i < width; i++)
	//{
	//	if (i == 0);   startPoint = Point2i(0, ptr[i]);
	//	if (i == width - 1) endPoint = Point2i(width - 1, ptr[i]);
	//	if (yvalue > ptr[i])
	//	{
	//		yvalue = ptr[i];
	//		vertexPoint = Point2i(i, yvalue);
	//	}
	//}

	////qDebug()<<"vertex Point is "<<vertexPoint.x<<","<<vertexPoint.y<<endl;
	//	//qDebug()<<"start Point is "<<startPoint.x<<","<<startPoint.y<<endl;
	//	//qDebug()<<"endPoint is"<<endPoint.x<<","<<endPoint.y<<endl;

	//double dist1 = sqrt((double)(vertexPoint.x - startPoint.x)*(vertexPoint.x - startPoint.x) + (vertexPoint.y - startPoint.y)*(vertexPoint.y - startPoint.y));
	//double dist2 = sqrt((double)(vertexPoint.x - endPoint.x)*(vertexPoint.x - endPoint.x) + (vertexPoint.y - endPoint.y)*(vertexPoint.y - endPoint.y));

	//Point2i featurePointA;
	//Point2i featurePointB;
	//double temp = 0;
	//double tempDist = 0;
	//double maxDist = 0;


	//Point2i vecA = vertexPoint - startPoint;
	//Point2i vecTemp;


	//for (int i = 0; i<vertexPoint.x; i++)
	//{
	//	if (ptr[i] == 0)
	//	{
	//		temp = 0;
	//		tempDist = 0;
	//		continue;
	//	}
	//	else
	//	{
	//		vecTemp = Point2i(i, ptr[i]) - startPoint;
	//		temp = abs(det(vecTemp, vecA));
	//		tempDist = temp / dist1;
	//	}

	//	if (maxDist < tempDist)
	//	{
	//		maxDist = tempDist;
	//		featurePointA = Point2i(i, ptr[i]);
	//	}
	//	
	//}

	//Point2i vecB = vertexPoint - endPoint;
	//temp = 0;
	//tempDist = 0;
	//maxDist = 0;
	//for (int i = vertexPoint.x; i < endPoint.x; i++)
	//{
	//	if (ptr[i] == 0)
	//	{
	//		temp = 0;
	//		tempDist = 0;
	//		continue;
	//	}
	//	else
	//	{
	//		vecTemp = Point2i(i, ptr[i]) - endPoint;
	//		temp = abs(det(vecTemp, vecB));
	//		tempDist = temp / dist2;
	//	}
	//	
	//	if (maxDist < tempDist)
	//	{
	//		maxDist = tempDist;
	//		featurePointB = Point2i(i, ptr[i]);
	//	}
	//}

	//featurePoints.push_back(featurePointA);
	//featurePoints.push_back(featurePointB);

	/*qDebug()<<"featurePointA"<<featurePoints[0].x<<","<<featurePoints[0].y<<endl;
		qDebug()<<"featurePointB"<<featurePoints[1].x<<","<<featurePoints[1].y<<endl;*/
	for (int i = 0; i < roiWidth - 1; i++)
	{
		slope = maxLoc[i] - maxLoc[i + 1];
		if (maxSlope < slope)
		{
			maxSlope = slope;
			maxSlopeLoc = i;
			featurePointA.x = i;
			featurePointA.y = maxLoc[i];
		}
		if (minSlope > slope)
		{
			minSlope = slope;
			minSlopeLoc = i;
			featurePointC.x = i;
			featurePointC.y = maxLoc[i];
		}
	}
	featurePointB.x = featurePointA.x + 1;
	featurePointB.y = maxLoc[featurePointB.x];
	featurePointD.x = featurePointC.x + 1;
	featurePointD.y = maxLoc[featurePointD.x];
	featurePoints.push_back(featurePointA);
	featurePoints.push_back(featurePointB);
	featurePoints.push_back(featurePointC);
	featurePoints.push_back(featurePointD);
	for (int i = 0; i <= maxSlopeLoc; i++)
	{
		t1 += i * i;
		t2 += i;
		t3 += maxLoc[i] * i;
		t4 += maxLoc[i];
	}
	b0 = (t1*t4 - t2 * t3) / (t1*maxSlopeLoc - t2 * t2);
	b1 = (t3*maxSlopeLoc - t2 * t4) / (t1*maxSlopeLoc - t2 * t2);
	for (int i = 0; i <= maxSlopeLoc; i++)
	{
		yi = b0 + b1 * i;
		point.x = ROIX +i;
		point.y = (int)(ROIY+yi);
		centerLine.push_back(point);
	}

	t1 = 0;
	t2 = 0;
	t3 = 0;
	t4 = 0;
	yi = 0;

	for (int i = (maxSlopeLoc + 1); i <= (minSlopeLoc); i++)
	{
		t1 += i * i;
		t2 += i;
		t3 += maxLoc[i] * i;
		t4 += maxLoc[i];
	}
	b0 = (t1*t4 - t2 * t3) / (t1*(minSlopeLoc - maxSlopeLoc) - t2 * t2);
	b1 = (t3*(minSlopeLoc - maxSlopeLoc) - t2 * t4) / (t1*(minSlopeLoc - maxSlopeLoc) - t2 * t2);


	for (int i = (maxSlopeLoc + 1); i <= (minSlopeLoc); i++)
	{
		yi = b0 + b1 * i;
		point.x = ROIX + i;
		point.y = (int)(ROIY + yi);
		centerLine.push_back(point);

	}

	t1 = 0; t2 = 0; t3 = 0; t4 = 0, yi = 0;
	for (int i = minSlopeLoc+1; i < roiWidth; i++)
	{
		t1 += i * i;
		t2 += i;
		t3 += maxLoc[i] * i;
		t4 += maxLoc[i];
	}
	b0 = (t1*t4 - t2 * t3) / (t1*(roiWidth - minSlopeLoc-1) - t2 * t2);
	b1 = (t3*(roiWidth - minSlopeLoc-1) - t2 * t4) / (t1*(roiWidth - minSlopeLoc-1) - t2 * t2);
	for (int i = minSlopeLoc+1; i < roiWidth; i++)
	{
		yi = b0 + b1 * i;
		point.x = ROIX + i;
		point.y = (int)(ROIY + yi);
		centerLine.push_back(point);
	}
	Mat leastSquare = Mat::zeros(Size(roiWidth, roiHeight), CV_8UC1);
	for (int i = 0; i < roiWidth; i++)
	{
		gray.at<uchar>(centerLine[i]) = 0;
	}

	/*time0 = ((double)getTickCount() - time0) / getTickFrequency();
	std::cout << "此方法运行时间为：" << time0 << "秒" << std::endl;*/
}


/*****************************************************************************
Function:void DIP::generateInformation()
Description:产生焊缝宽度、焊缝深度、焊缝偏差信息
Input:
Output:
Return:
Others:
*****************************************************************************/
void DIP::generateInformation()
{
	double centerx;
	centerx = (featurePoints[1].x + featurePoints[2].x) / 2;
	offset = roiWidth / 2 - centerx;
	weldSeamWidth = (featurePoints[2].x + featurePoints[3].x) / 2 - (featurePoints[0].x + featurePoints[1].x) / 2;
	weldSeamDepth = (featurePoints[0].y + featurePoints[3].y) / 2 - (featurePoints[1].y + featurePoints[2].y) / 2;
	std::cout << "焊缝宽度：" << weldSeamWidth << std::endl;
	std::cout << "焊缝深度：" << weldSeamDepth << std::endl;
	std::cout << "焊缝偏差：" << offset << std::endl;
}


/*****************************************************************************
Function:void DIP::generateOutImage()
Description:产生Out图，此图片含有拟合出来的直线和焊缝特征点
Input:
Output:
Return:
Others:
*****************************************************************************/



void DIP::generateOutImage()
{
	//用小十字架标出特征点位置
	int horizontalShift = (width - roiWidth) / 2;
	Point2i featurePointA = Point2i(ROIX + featurePoints[0].x, featurePoints[0].y + ROIY);
	Point2i featurePointB = Point2i(ROIX + featurePoints[1].x, featurePoints[1].y + ROIY);
	Point2i featurePointC = Point2i(ROIX + featurePoints[2].x, featurePoints[2].y + ROIY);
	Point2i featurePointD = Point2i(ROIX + featurePoints[3].x, featurePoints[3].y + ROIY);
	drawAsterisk(gray, featurePointA, Scalar(0,0, 0), 2);
	drawAsterisk(gray, featurePointB, Scalar(0, 0, 0), 2);
	drawAsterisk(gray,featurePointC, Scalar(0, 0, 0), 2);
	drawAsterisk(gray, featurePointD, Scalar(0, 0, 0), 2);


	////用小十字架标出焊缝中心位置
	//Point2i absWeldSeamCenterPoint = 0.5*(leftFeaturePoint + rightFeaturePoint);
	//drawAsterisk(image, absWeldSeamCenterPoint, Scalar(255, 0, 0), 2);

	//用蓝色线框标出结构光条纹图像ROI的位置
	rectangle(gray, Point2i(ROIX, ROIY), Point2i(ROIX + roiWidth, ROIY + roiHeight), Scalar(255, 0, 0), 2, 8);

	////用蓝色线画出中心线的位置
	//startPoint.x = roiX;
	//startPoint.y = roiY+maxLoc[0];
	//endPoint.x = roiX + roiWidth - 1;
	//endPoint.y = roiY + maxLoc[roiWidth - 1];
	//line(image, startPoint, featurePointA, Scalar(255, 0, 0), 2, 8, 0);
	//line(image, featurePointA, featurePointB, Scalar(255, 0, 0), 2, 8, 0);
	//line(image, featurePointB, featurePointC, Scalar(255, 0, 0), 2, 8, 0);
	//line(image, featurePointC, featurePointD, Scalar(255, 0, 0), 2, 8, 0);
	//line(image, featurePointD, endPoint, Scalar(255, 0, 0), 2, 8, 0);
	////用蓝色实线标出焊枪中心的位置
	//line(image, Point2i(roiX + roiWidth / 2, roiY), Point2i(roiX + roiWidth / 2, roiY + roiHeight), Scalar(255, 0, 0), 2, 8);

	////标示“ROI”字样和焊缝偏差信息
	//String roiTitle = "ROI";
	//putText(image, roiTitle, Point2i(roiX, roiY - 5), FONT_HERSHEY_SIMPLEX, 1.5, Scalar(0, 0, 255), 2, 8);

	//String offsetInfo = "Offset" + std::to_string(long double(offset)) + "pixels";
	//putText(image, offsetInfo, Point2i(20, height - 40), FONT_HERSHEY_SIMPLEX, 1.5, Scalar(0, 0, 255), 2, 8);

	////标示焊缝宽度
	//String weldSeamWidthInfo = "WeldSeamWidth" + std::to_string(long double(weldSeamWidth)) + "pixels";
	//putText(image, weldSeamWidthInfo, Point2i(20, height - 80), FONT_HERSHEY_SIMPLEX, 1.5, Scalar(0, 0, 255), 2, 8);

	outputImage = gray;
	this->out = DIP::mat2QImage(gray);

}  


/*****************************************************************************
Function:mat2QImage
Description:Mat转QImage
Input:
Output:
Return:
Others:
*****************************************************************************/
QImage DIP::mat2QImage(Mat &in)
{
    if(in.type()==CV_8UC1)
    {
        QImage img((const unsigned char *)(in.data),in.cols,in.rows,in.step,QImage::Format_Grayscale8);
        return img;
    }
    else if(in.type()==CV_8UC3)
    {
        const uchar *pSrc = (const uchar *)in.data;
        QImage out(pSrc,in.cols,in.rows,in.step,QImage::Format_RGB888);
        return out.rgbSwapped();
    }
}

/*****************************************************************************
Function:det
Description:求两个向量的行列式
Input:
Output:
Return:
Others:
*****************************************************************************/
double DIP::det(Point2i a, Point2i b)
{
	double res;
	res = a.x*b.y - a.y*b.x;
	return res;
}

/*****************************************************************************
Function:drawAsterisk
Description:画特征点的星号
Input:
Output:
Return:
Others:
*****************************************************************************/
void DIP::drawAsterisk(Mat& image, Point pt, Scalar& color, int thickness)
{
	line(image, Point2i(pt.x - 3, pt.y), Point2i(pt.x + 3, pt.y), color, thickness);
	line(image, Point2i(pt.x, pt.y-5), Point2i(pt.x, pt.y+5), color, thickness);
}