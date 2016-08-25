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

class MaxAbsScalerObject {
public:
	int label;
	std::vector<double> values;
	double e_max;
	double e_min;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}
}

class MaxAbsScaler{
private:
	std::string inputCol;
	std::string outputCol;
	double min;
	double max;

	Husky::ObjList<MaxAbsScalerObject>& inputList;
	
public:
	MaxAbsScaler(){
		min = -1.0;
		max = 1.0;
		inputCol = "inputVector";
		outputCol = "outputVector";
	};

	~MaxAbsScaler();

	Husky::ObjList<MaxAbsScalerObject>& getInputList() {
		return inputList;
	};

	std::string getInputCol(){
		return inputCol;
	};

	std::string getOutputCol() {
		return outputCol;
	};

	double getMin() {
		return min;
	};

	double getMax() {
		return max;
	};

	void setInputCol(std::string input){
		inputCol = input;
	};

	void setOutputCol(std::string output){
		outputCol = output;
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
				MaxAbsScalerObject tempMaxAbsScalerObject;
				if (marker % 2 == 0) {
					boost::char_separator<char> sep2(",");
					boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);
					for (auto &c : tok2){
						tempMaxAbsScalerObject.values.push_back(std::stod(c));
					}
					worker.add_object(inputList, tempMaxAbsScalerObject);
				}
				else{
					tempMaxAbsScalerObject.label = std::stoi(b);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};

	void transform(){

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		worker.list_execute(inputList, [&](MaxAbsScalerObject &tmpObj){
			int len = tmpObj.values.size();

			tmpObj.e_min = INT_MAX;
			tmpObj.e_max = 0;

			for (int i = 0; i < len; i++) {
				if (tmpObj.values[i] > tmpObj.e_max) {
					tmpObj.e_max = tmpObj.values[i];
				}

				if (tmpObj.values[i] < tmpObj.e_min) {
					tmpObj.e_min = tmpObj.values[i];
				}
			}

			if (tmpObj.e_max == tmpObj.e_min){
				for (auto &x : tmpObj.values){
					x = 0.5 * (max + min);
				}
			}
			else {
				for (auto &x : tmpObj.values){
					x = ((x - e_min) / (e_max - e_min)) * (max - min) + min;
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};
}