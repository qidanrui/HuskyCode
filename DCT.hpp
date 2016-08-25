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

class DCTObject {
public:
	int label;
	std::vector<double> features;
	std::vector<double> dctFeatures;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}
}

class DCT {
private:
	std::string inputCol;
	std::string outputCol;

	Husky::ObjList<DCTObject>& inputList;
	
public:
	DCT(){
		inputCol = "inputTokens";
		outputCol = "nGram";
	};

	DCT(int para, std::string inputColumn, std::string outputColumn){
		inputCol = inputColumn;
		outputCol = outputColumn;
	};
	~DCT();

	Husky::ObjList<DCTObject>& getInputList() {
		return inputList;
	};

	std::string getInputCol(){
		return inputCol;
	};

	std::string getOutputCol() {
		return outputCol;
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
				DCTObject tempDCTObject;
				if (marker % 2 == 0) {
					boost::char_separator<char> sep2(",");
					boost::tokenizer<boost::char_separator<char>> tok2(b, sep2);
					for (auto &c : tok2){
						tempDCTObject.features.push_back(c);
					}
					worker.add_object(inputList, tempDCTObject);
				}
				else{
					tempDCTObject.label = std::stoi(b);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};

	void transform(){

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		worker.list_execute(inputList, [&](DCTObject &tmpObj){
			len = tmpObj.features.size();
			for (int i = 0; i < len; i++) {
				double res = 0;
				for (int j = 0; j < len; j++){
					res += tmpObj.features[j] * cos((pi * (j + 0.5) * i) / len);
				}
				tmpObj.dctFeatures.push_back(res);
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};
}