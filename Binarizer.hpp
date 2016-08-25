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

class BinarizerObject {
public:
	int label;
	std::vector<double> inputValues;
	std::vector<double> features;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}
}

class Binarizer {
private:
	std::string inputCol;
	std::string outputCol;
	double threshold;

	Husky::ObjList<NGramObject>& inputList;
	
public:
	Binarizer(){
		threshold = 1;
		inputCol = "inputTokens";
		outputCol = "nGram";
	};

	Binarizer(double para, std::string inputColumn, std::string outputColumn){
		threshold = para;
		inputCol = inputColumn;
		outputCol = outputColumn;
	};
	~Binarizer();

	Husky::ObjList<NGramObject>& getInputList() {
		return inputList;
	};

	std::string getInputCol(){
		return inputCol;
	};

	std::string getOutputCol() {
		return outputCol;
	};

	double getThreshold() {
		return threshold;
	}

	void setInputCol(std::string input){
		inputCol = input;
	};

	void setOutputCol(std::string output){
		outputCol = output;
	};

	void setParams(double para){
		threshold = para;
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
				BinarizerObject tempBinarizerObject;
				if (marker % 2 == 0) {
					boost::char_separator<char> sep2(",");
					boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);
					for (auto &c : tok2){
						tempBinarizerObject.inputValues.push_back(c);
					}
					worker.add_object(inputList, tempBinarizerObject);
				}
				else{
					tempBinarizerObject.label = std::stod(b);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};

	void transform(){

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		worker.list_execute(inputList, [&](BinarizerObject &tmpObj){
			for (auto &x : tmpObj.inputValues) {
				if (x >= threshold) {
					tmpObj.features.push_back(1.0);
				}
				else {
					tmpObj.features.push_back(0.0);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};
}

void runWorker() {
	Husky::Basewoker & worker = Husky::Context::get_worker<Husky::Basewoker>();
	auto & binarizerList = worker.create_list<BinarizerObject>("binarizer_list");

	Binarizer binarizer(binarizerList);

	binarizer.loadFile(worker, "input_binarizer");
	binarizer.transform(worker);
}

int main(int argc, char **argv) {
	Husky::run_job(runWorker, argv[1]);
	return 0;
}