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

class StandardScalerObject {
public:
	int label;
	std::vector<double> values;
	double mean;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}
}

class StandardScaler{
private:
	std::string inputCol;
	std::string outputCol;
	bool withStd;
	bool withMean;

	Husky::ObjList<StandardScalerObject>& inputList;
	
public:
	StandardScaler(){
		withStd = true;
		withMean = false;
		inputCol = "inputVector";
		outputCol = "outputVector";
	};

	StandarScaler(bool para1, bool para2, std::string inputColumn, std::string outputColumn){
		withStd = para1;
		withMean = para2;
		inputCol = inputColumn;
		outputCol = outputColumn;
	};

	~StandardScaler();

	Husky::ObjList<StandardScalerObject>& getInputList() {
		return inputList;
	};

	std::string getInputCol(){
		return inputCol;
	};

	std::string getOutputCol() {
		return outputCol;
	};

	bool getWithStd() {
		return withStd;
	};

	bool getWithMean() {
		return withMean;
	};

	void setInputCol(std::string input){
		inputCol = input;
	};

	void setOutputCol(std::string output){
		outputCol = output;
	};

	void setWithStd(bool para){
		withStd = para;
		withMean = !para;
	};

	void setWithMean(bool para){
		withStd = !para;
		withMean = para;
	}

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
				StandardScalerObject tempStandardScalerObject;
				if (marker % 2 == 0) {
					boost::char_separator<char> sep2(",");
					boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);
					for (auto &c : tok2){
						tempStandardScalerObject.values.push_back(std::stod(c));
					}
					worker.add_object(inputList, tempStandardScalerObject);
				}
				else{
					tempStandardScalerObject.label = std::stoi(b);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};

	void transform(){

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		worker.list_execute(inputList, [&](StandardScalerObject &tmpObj){
			int len = tmpObj.values.size();

			tmpObj.mean = 0;

			for (int i = 0; i < len; i++) {
				tmpObj.mean += tmpObj.values[i];
			}

			tmpObj.mean /= len;

			if (withMean) {

				for (auto &x : tmpObj.values){
					x -= tmpObj.mean;
				}
			}
			else {
				x *= sqrt(len / tmpObj.mean);
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};
}