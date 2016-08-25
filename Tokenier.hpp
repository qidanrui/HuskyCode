#include <algorithm>
#include <cfloat>
#include <string>
#include <vector>
#include <functional>

#include "boost/tokenizer.hpp"
#include "core/context.hpp"
#include "core/engine.hpp"
#include "lib/dcaaggregator.hpp"
#include "lib/topk.hpp"
#include "gperftools/profiler.h"
#include "modules/inputformat/lineinputformat.hpp"

class TokenierObject {
public:
	int label;
	std::string sentence;
	std::vector<std::string> words;
	typedef int KeyT;

	virtual KeyT const & id() const {
		return label;
	}
}

class Tokenier : public Husky::BaseWorker {
private:
	std::string inputCol;
	std::string outputCol;

	Husky::ObjList<TokenierObject>& inputList;
	
public:
	Tokenier();
	Tokenier(std::string inputColumn, std::string outputColumn);
	~Tokenier();

	Husky::ObjList<TokenierObject>& getInputList() {
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
				TokenierObject tempTokenier;
				if (marker % 2 == 0) {
					TokenierObject.sentence = b;
					worker.add_object(inputList, tempTokenier);
				}
				else{
					tempTokenier.label = std::stoi(b);
				}
			}
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};

	void transform(){
		regex("\\W");
	};
	
	void regex(std::string pattern){
		auto worker = Husky::Context::get_worker<Husky::Basewoker>();

		worker.list_execute(inputList, [&](TokenierObject &tObj){
			boost::regex re(pattern);
			boost::sregex_token_iterator i(tObj.sentence.begin(), tObj.sentence.end(), re, -1);
			boost::sregex_token_iterator j;

			while(i != j)
			{
				tObj.words.push_back(*i++); 
			}
 
		});

		Husky::Context::free_worker<Husky::BaseWorker>();
	};
}

void runWorker(){
	Tokenier & worker = Husky::Context::get_context().get_worker<Tokenier>();
	worker.loadFile("input");
	worker.regex(",");
	worker.transform();
    Husky::Context::get_context().free_worker<Tokenier>();
}

int main(int argc, char** argv) {
    Husky::run_job(runWorker, argv[1]);
    return 0;
}

