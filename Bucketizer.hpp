#include <algorithm>
#include <cfloat>
#include <string>
#include <vector>
#include <functional>

#include "core/context.hpp"
#include "core/engine.hpp"
#include "lib/dcaaggregator.hpp"
#include "lib/topk.hpp"
#include "gperftools/profiler.h"
#include "modules/inputformat/lineinputformat.hpp"

class BucketizerObject {
public:
	int label;
	std::vector<double> inputValues;
	std::unordered_map<std::pair<double, double>, std::vector<double>> bucket;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}
};

class Bucketizer {
private:
	std::string inputCol;
	std::string outputCol;
	std::vector<double> bucketizer;

	Husky::ObjList<BucketizerObject>& inputList;
	
public:
	Bucketizer(){
		inputCol = "inputTokens";
		outputCol = "Bucketizer";
	};

	Bucketizer(std::vector<double> para, std::string inputColumn, std::string outputColumn){
		bucketizer.clear();
		for (auto &x ：para) {
			bucketizer.push_back(x);
		}
		inputCol = inputColumn;
		outputCol = outputColumn;
	};
	~Bucketizer();

	Husky::ObjList<BucketizerObject>& getInputList() {
		return inputList;
	};

	std::string getInputCol(){
		return inputCol;
	};

	std::string getOutputCol() {
		return outputCol;
	};

	std::vector<double> getBucketizer(){
		return bucketizer;
	};

	void setInputCol(std::string input){
		inputCol = input;
	};

	void setOutputCol(std::string output){
		outputCol = output;
	};

	void setBucketizer(std::vector<double> para){
		bucketizer.clear();
		for (auto &x : para){
			bucketizer.push_back(x);
		}
	};

	bool sortFunc(double v1, double v2)//注意：本函数的参数的类型一定要与vector中元素的类型一致  
	{  
		return v1 < v2;//升序排列  
	};

	void loadFile(std::string fileNamePara) {
		std::string input_col = inputCol;
		std::stirng output_col = outputCol;

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		Husky::HDFSLineInputFormat infmt;
		infmt.set_input(Husky::Context::get_params(fileNamePara));

		int marker = 0;

		worker.load(infmt, [&](boost::string_ref chunk) {
			if (chunk.empty()) return;
			boost::char_separator<char> sep(" \t");
			boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);

			for (auto &b : tok) {
				marker++;
				BucketizerObject tempBucketizerObject;
				if (marker % 2 == 0) {
					boost::char_separator<char> sep2(",");
					boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);
					for (auto &c : tok2){
						tempBucketizerObject.inputValues.push_back(std::stod(c));
					}
					std::sort(tmpObj.inputValues.begin(), tmpObj.inputValues.end(), sortFunc);  
					worker.add_object(inputList, tempBucketizerObject);
				}
				else{
					tempBucketizerObject.label = std::stoi(b);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};

	void transform(){

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		worker.list_execute(inputList, [&](BucketizerObject &tmpObj){
			int bucketizerLength = bucketizer.size();
			int inputValuesLength = tmpObj.inputValues.size();
			int i = 0;
			int j = 0;
			while (i < bucketizerLength - 1){
				while (j < inputValuesLength){
					if ((tmpObj.inputValues[j] >= bucketizer[i]) && (tmpObj.inputValues[j] < bucketizer[i + 1])){
						tmpObj.bucket[std::make_pair(bucket[i], bucket[i + 1])].push_back(tmpObj.inputValues[j]);
						j++;
					}
					else{
						i++;
						break;
					}
				}
			}	
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};
}