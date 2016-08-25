//输入文件 已经分好的词
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

class SWRObject {
public:
	std::string word;
	
	typedef string KeyT;

	virtual KeyT const & id() const {
		return word;
	}
}

class StopWordsRemover {
private:
	bool caseSensitive;
	std::set<std::string> stopWords;
	std::vector<std::string> res;

	Husky::ObjList<SWRObject>& inputList;
	
public:
	StopWordsRemover() : inputList(inputList){
		caseSensitive = false;
	};

	virtual ~StopWordsRemover();

	Husky::ObjList<TokenierObject>& getInputList() {
		return inputList;
	};	

	bool getCaseSensitive(){
		return caseSensitive;
	}

	void setCaseSensitive(bool sensi){
		caseSensitive = sensi;
	}

	void loadDefaultStopWords(Husky::BaseWorker &worker, std::string language) {

		Husky::HDFSLineInputFormat infmt;
		infmt.set_input(Husky::Context::get_params("DefaultStopWords-" + language));
		int marker = 0;

		worker.load(infmt, [&](boost::string_ref chunk) {
			if (chunk.empty()) return;
			boost::char_separator<char> sep(" \t");
			boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);
			for (&b : tok){
				stopWords.insert(b);
			}
			/*SWRObject tempStopWordsRemover;
			tempStopWordsRemover.word = chunk.data();*/
		});
	}

	void addStopWords(Husky::BaseWorker &worker, std::string inputStopWords){
		Husky::HDFSLineInputFormat infmt;
		infmt.set_input(Husky::Context::get_params(inputStopWords));

		worker.load(infmt, [&](boost::string_ref chunk) {
			if (chunk.empty()) return;
			boost::char_separator<char> sep(" \t");
			boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);
			for (&b : tok){
				stopWords.insert(b);
			}
			/*SWRObject tempStopWordsRemover;
			tempStopWordsRemover.word = chunk.data();*/
		});
	}

	void setStopWords(Husky::BaseWorker &worker, std::string inputStopWords){
		stopWords.clear();

		Husky::HDFSLineInputFormat infmt;
		infmt.set_input(Husky::Context::get_params(inputStopWords));

		worker.load(infmt, [&](boost::string_ref chunk) {
			if (chunk.empty()) return;
			boost::char_separator<char> sep(" \t");
			boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);
			for (&b : tok){
				stopWords.insert(b);
			}
			/*SWRObject tempStopWordsRemover;
			tempStopWordsRemover.word = chunk.data();*/
		});
	}

	void initStopWords(Husky::BaseWorker &worker, std::string language){
		stopWords.clear();
		loadDefaultStopWords(worker, language);		
	}

	void loadInputFile(HUsky::BaseWorker &worker, std::string fileName) {

		Husky::HDFSLineInputFormat infmt;
		infmt.set_input(Husky::Context::get_params(fileNamePara));
		int marker = 0;

		worker.load(infmt, [&](boost::string_ref chunk) {
			if (chunk.empty()) return;
			boost::char_separator<char> sep(" \t");
			boost::tokenizer<boost::char_separator<char>> tok(chunk, sep);
			for (&b : tok){
				SWRObject tempStopWordsRemover;
				tempStopWordsRemover.word = b;
				worker.add_object(inputList, tempStopWordsRemover);
			}
		});
	}

	void transform(Husky::BaseWorker &worker, bool caseSensi){
		setCaseSensitive(caseSensi);
		if (caseSensitive){
			for (auto &x : stopWords) {
				x = strlwr(x);
			}

			worker.list_execute(inputList, [&](SWRObject &tmpObj){
				tmpObj.word = strlwr(tmpObj.word);
			});
		}

		worker.list_execute(inputList, [&](SWRObject &tmpObj){
			if (stopWords.find(tmpObj.word) == stopWords.end()) {
				continue;
			}
			else {
				res.push_back(tmpObj.word);
			}
		});
	}
};

void runWorker(){
	Husky::BaseWorker & worker = Husky::Context::get_worker<Husky::BaseWorker>();
	auto & stopWordsRemoverList = worker.create_list<SWRObject>("stopWordsRemover_list");

	StopWordsRemover remover;
	remover.addStopWords(worker, "inputStopWords");
	remover.loadInputFile(worker, "inputFile");
	remover.transform(worker, true);
    Husky::Context::get_context().free_worker<Husky::BaseWorker()>();
}

int main(int argc, char** argv) {
    Husky::run_job(runWorker, argv[1]);
    return 0;
}

