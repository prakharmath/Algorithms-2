//Write a C/C++ program to generate n random points, find the closest pair in O(nlogn) time, and draw discs of same radius centred at these points so that the discs at the closest pair just touch each other. The output should be an SVG file.


#include<bits/stdc++.h>
using namespace std;
float edis(pair<int, int> a, pair<int, int> b)
{
	int xval = (abs(a.first - b.first) * abs(a.first - b.first));
	int yval = (abs(a.second - b.second) * abs(a.second - b.second));
	return sqrt((xval + yval));
}
bool sortsecond(pair<int, int> a, pair<int, int> b)
{
	if (a.second < b.second)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
float between(vector<pair<int, int>> b, float min_val)
{
	float answer = min_val;
	sort(b.begin(), b.end(), sortsecond);
	for (auto itr = b.begin(); itr != b.end(); itr++)
	{
		// This inside loop runs at most 6 times
		for (auto itr2 = itr + 1; ((itr2 != b.end()) && ((itr2->second - itr->second) <= answer)); itr2++)
		{
			auto one_point = make_pair(itr2->first, itr2->second);
			auto other_point = make_pair(itr->first, itr->second);
			float distance_between = edis(one_point, other_point);
			if (distance_between < answer)
			{
				answer = distance_between;
			}
		}
	}
	// printf("From combining the answer is %f\n",answer);
	return answer;
}
float closestpair(vector<pair<int, int>> v, int l, int r)
{
	if (l == r)
	{
		return INT_MAX;
	}
	int mid = (l + r) / 2;
	float left_side = closestpair(v, l, mid);
	float right_side = closestpair(v, mid + 1, r);
	float min_val = min(left_side, right_side);
	vector<pair<int, int>> b;
	for (int i = l; i <= r; i++)
	{
		// taking a rectangle of width 2*min_val
		if (abs(v[i].first - v[mid].first) <= min_val)
		{
			b.push_back(v[i]);
		}
	}
	float value = between(b, min_val);
	min_val = min(value, min_val);
	// printf("For %d %d value is %f\n",l,r,min_val);
	return min_val;
}
int main()
{
	int n, i;
	printf("Write the number of points\n");
	scanf("%d", &n);
	vector<pair <int, int>> v;
	printf("Points will be randomely generated ranged from 100 to 299\n The canvas width and lenght both are 1000\n");
	int x, y;
	for (i = 0; i < n; i++)
	{
		x = 100 + rand() % 300;
		y = 100 + rand() % 300;
		v.push_back(make_pair(x, y));
	}
	sort(v.begin(), v.end());
	printf("The points in sorted order ( according to x coordinate ) are:\n");
	for (int i = 0; i < n; i++)
	{
		printf("%d %d\n", v[i].first, v[i].second);
	}
	float answer = closestpair(v, 0, n - 1);
	printf("The closest distance is %f\n", answer);
	float ranswer = answer / 2; // radius is half the distance
	// Now start making the svg file
	ofstream outdata;
	// svg start
	string start = "<!DOCTYPE html><html><body><svg width=\"1000\" height=\"1000\">";
	// circle svg code
	string circle_start = "<circle ";
	string circle_attributes = "fill=\"yellow\" fill-opacity=\"0.4\" ";
	string circle_end = "/>";
	string radius = "r=";
	radius.append("\"");
	radius.append(to_string(ranswer));
	radius.append("\" ");
	// text svg code
	string text_start = "<text ";
	string text_attributes = "fill=\"blue\">";
	string text_end = "</text>";
	// svg end code
	string endline = "</svg></body></html>";
	// open output file
	outdata.open("image.svg");
	if (!outdata)
	{
		printf("Error occurred while opening file, aborting!");
		return 0;
	}
	outdata << start;
	for (int j = 0; j < n; j++)
	{
		// format the circle string
		string cxcoordinate = "cx=\"";
		cxcoordinate.append(to_string(v[j].first));
		cxcoordinate.append("\" ");
		string cycoordinate = "cy=\"";
		cycoordinate.append(to_string(v[j].second));
		cycoordinate.append("\" ");
		string circle = "";
		circle.append(circle_start);
		circle.append(cxcoordinate);
		circle.append(cycoordinate);
		circle.append(radius);
		circle.append(circle_attributes);
		circle.append(circle_end);
		// format the text string
		string txcoordinate = "x=\"";
		txcoordinate.append(to_string(v[j].first));
		txcoordinate.append("\" ");
		string tycoordinate = "y=\"";
		tycoordinate.append(to_string(v[j].second));
		tycoordinate.append("\" ");
		string text = "";
		text.append(text_start);
		text.append(txcoordinate);
		text.append(tycoordinate);
		text.append(text_attributes);
		text.append(to_string(j + 1));
		text.append(text_end);
		outdata << circle << text;
	}
	outdata << endline;
	outdata.close();
	return 0;
}
