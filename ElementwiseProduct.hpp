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

class ElementwiseProductObject {
public:
	int label;
	std::vector<double> inputValues;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}
};

class ElementwiseProduct {
private:
	std::string inputCol;
	std::string outputCol;
	std::vector<double> transformVector;

	Husky::ObjList<ElementwiseProductObject>& inputList;
	
public:
	ElementwiseProduct(){
		inputCol = "inputTokens";
		outputCol = "ElementwiseProduct";
	};

	Bucketizer(std::vector<double> para, std::string inputColumn, std::string outputColumn){
		transformVector.clear();
		for (auto &x ：para) {
			transformVector.push_back(x);
		}
		inputCol = inputColumn;
		outputCol = outputColumn;
	};

	~ElementwiseProduct();

	Husky::ObjList<ElementwiseProductObject>& getInputList() {
		return inputList;
	};

	std::string getInputCol(){
		return inputCol;
	};

	std::string getOutputCol() {
		return outputCol;
	};

	std::vector<double> getTransformVector(){
		return transformVector;
	};

	void setInputCol(std::string input){
		inputCol = input;
	};

	void setOutputCol(std::string output){
		outputCol = output;
	};

	void setTransformVector(std::vector<double> para){
		transformVector.clear();
		for (auto &x : para){
			transformVector.push_back(x);
		}
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
				ElementwiseProductObject tempElementwiseProductObject;
				if (marker % 2 == 0) {
					boost::char_separator<char> sep2(",");
					boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);
					for (auto &c : tok2){
						tempElementwiseProductObject.inputValues.push_back(std::stod(c));
					}
					std::sort(tmpObj.inputValues.begin(), tmpObj.inputValues.end(), sortFunc);  
					worker.add_object(inputList, tempElementwiseProductObject);
				}
				else{
					tempElementwiseProductObject.label = std::stoi(b);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};

	void transform(){

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		worker.list_execute(inputList, [&](ElementwiseProductObject &tmpObj){
			int len = tmpObj.inputValues.size();

			for (int i = 0; i < len; i++) {
				tmpObj.inputValues[i] *= transformVector[i];
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};
}