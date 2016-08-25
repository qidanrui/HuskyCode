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

class NGramObject {
public:
	int label;
	std::vector<std::string> inputTokens;
	std::vector<std::string> nGrams;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}
}

class nGram {
private:
	std::string inputCol;
	std::string outputCol;
	int n;

	Husky::ObjList<NGramObject>& inputList;
	
public:
	StopWordsRemover(){
		n = 1;
		inputCol = "inputTokens";
		outputCol = "nGram";
	};

	StopWordsRemover(int para, std::string inputColumn, std::string outputColumn){
		n = para;
		inputCol = inputColumn;
		outputCol = outputColumn;
	};
	~StopWordsRemover();

	Husky::ObjList<NGramObject>& getInputList() {
		return inputList;
	};

	int getN() {
		return n;
	}

	void setParams(int para){
		n = para;
	}

	void loadFile(std::string fileNamePara) {

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		Husky::HDFSLineInputFormat infmt;
		infmt.set_input(Husky::Context::get_params(fileNamePara));

		int marker = 0;

		worker.load(infmt, [&](boost::string_ref chunk) {
			if (chunk.empty()) return;
			boost::char_separator<char> sep(",");
			boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);
			for (auto &c : tok) {
				NGramObject tempNGramObject;
				tempNGramObject.inputTokens.push_back(c);
				worker.add_object(inputList, tempNGramObject); 
			}
		});
	};

	void transform(){

		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		worker.list_execute(inputList, [&](NGramObject &tmpObj){
			len = tmpObj.inputTokens.size();
			for (int i = 0; i < len - n + 1; i++) {
				std::string tmpNGramStr = "";
				for (int j = 0; j < n; j++){
					tmpNGramStr += tmpObj.inputTokens[i + j];
				}
				tmpObj.push_back(tmpNGramStr);
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};
}