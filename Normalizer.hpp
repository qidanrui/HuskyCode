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

class NormalizerObject {
public:
	int label;
	std::vector<double> values;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}
}

class Normalizer{
private:
	std::string inputCol;
	std::string outputCol;
	double p;

	Husky::ObjList<NormalizerObject>& inputList;
	
public:
	Normalizer(){
		p = 2.0;
		inputCol = "inputVector";
		outputCol = "outputVector";
	};

	Normalizer(double para, std::string inputColumn, std::string outputColumn){
		p = para;
		inputCol = inputColumn;
		outputCol = outputColumn;
	};
	~Normalizer();

	Husky::ObjList<NormalizerObject>& getInputList() {
		return inputList;
	};

	std::string getInputCol(){
		return inputCol;
	};

	std::string getOutputCol() {
		return outputCol;
	};

	double getP() {
		return p;
	}

	void setInputCol(std::string input){
		inputCol = input;
	};

	void setOutputCol(std::string output){
		outputCol = output;
	};

	void setParams(double para){
		p = para;
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
				NormalizerObject tempNormalizerObject;
				if (marker % 2 == 0) {
					boost::char_separator<char> sep2(",");
					boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);
					for (auto &c : tok2){
						tempNGramObject.values.push_back(std::stod(c));
					}
					worker.add_object(inputList, tempNormalizerObject);
				}
				else{
					tempNormalizerObject.label = std::stoi(b);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};

	void transform(){

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		worker.list_execute(inputList, [&](NormalizerObject &tmpObj){
			int len = tmpObj.inputValues.size();
			double dis = 0;
			for (int i = 0; i < len; i++) {
				dis += std::pow(tmpObj.values[i], p);
			}

			dis = pow(dis, 1.0 / p);

			for (auto &x : tmpObj.values){
				x = x / dis;
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};
}