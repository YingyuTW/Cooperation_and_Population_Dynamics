#include <iostream>
#include <math.h>
#include <vector>
#include <random>
#include <chrono>
using namespace std;

//Function to compute standard error
double StandardError(vector<double> data, double mean);

//To generate random number
class RandGenerator{
public:
	RandGenerator();
    double	RUnif();            //random generation for the continuous uniform distribution on the interval [0,1]
	double	Rand2();            //discrete  [0,1]
    bool	RBern(double p);    //random generation for the Bernoulli distribution
	int		SampleInt(int start, int end);  //random generation for the discrete uniform distribution
private:
	unsigned seed;
	default_random_engine RAND_GENERATOR;
};

//To create individuals
class Individual{
public:
    Individual(int PatchID, double CoopDeg, int Age);
	void	Aging(){ _Age = _Age + 1; }
    void	WriteReproRate(int    ReproRate){ _ReproRate = ReproRate; }
	void	WriteSurvRate (double SurvRate ){ _SurvRate  = SurvRate;  }
    double	GetCoopDeg()    { return _CoopDeg;     }
	double	GetSurvRate()   { return _SurvRate;    }
    int		GetReproRate()  { return _ReproRate;   }
	int		GetAge()        { return _Age;         }

private:
	int		_Age;
    int     _PatchID;
    int     _ReproRate;
	double	_SurvRate;
    double  _CoopDeg;
	RandGenerator RG;
};

//To create patches (social groups)
class Patch{
public:
	Patch();
    void	AddIndividual(int PatchID, double CoopDeg, int Age);
	void	UpdateAge();
    void	UpdateCoopNum();
    void	UpdateCoopDegSum();
    void	UpdateResource (double InitialResource, double CoopEfficiency, double MaxIncreRate);
    void	UpdateReproRate(double MetabolicConsump, double HalfConst, double ReproMax, double CostRate, bool End);
	void	UpdateSurvRate (double SurvRateUL, double AgeStandard, double CostRate);
	void	UpdateAlive();
	int		GetSize()		   		  { return _Individual.size();  }
    int		GetCoopNum()   		      { return _CoopNum;            }
	double	GetResource()	   		  { return _Resource;           }
    double	GetOffspringNum()  		  { return _OffspringNum;	    }
    double	GetCoopDegSum()		      { return _CoopDegSum;	        }
	//void	EraseIndividual    (int i){	_Individual.erase(_Individual.begin()+i);	}
	vector<double> GetOffspring();
    vector<double> GetEveryIndivRepro(){ return _EveryIndivRepro;   }
	vector<double> GetEveryIndivDeg()  { return _EveryIndivDeg;		}

private:
    int     _OffspringNum;
    int     _CoopNum;
    double  _CoopDegSum;
	double	_Resource;
    vector<double> _EveryIndivRepro;
	vector<double> _EveryIndivDeg;
	vector<Individual> _Individual;
	RandGenerator RG;
};

//To manage processes in population level
class PopProcess{
public:
	PopProcess(double CoopEfficiency, double EnvironResDens, double CostRate, double CoopProp, int Span, int CutPoint, double FlucAmp, double FlucPeriod, bool Social, bool Divided);
	~PopProcess();
    void    ComputePatchReproData();
    void	ComputePopCoopData();
    void    ComputeTimeSeriesData(int index, double FlucRes);
	void	OffspringDisperse();
    void	SystemUpdate(int Years);
	void	ErrorHandling1();	//when pop. size = 0
	void	ErrorHandling2();	//when pop. size is too large
    int     GetLastNA()                 { return _LastNA;           }
	int		GetPopSize()				{ return _PopSize;	        }
    int		GetTotalCoopNum()   		{ return _TotalCoopNum;	    }
    double	GetCoopProp()   			{ return _CoopProp;		    }
    double	GetAverCoopDeg()    		{ return _AverCoopDeg;	    }
	vector<double> GetTS()				{ return _TSData;			}
    vector<double> GetPatchReproData()	{ return _PatchReproData;	}
    vector<double> GetIndivReproData()	{ return _IndivReproData;	}
	vector<double> GetIndivDegData()	{ return _IndivDegData;		}

private:
	//Population state
    int     _PopSize;
	int		_PatchNum;
    int		_PopSizeUL;
    int		_TotalCoopNum;
    bool    _Divided;
    double	_CoopProp;
    double	_AverCoopDeg;
    double  _PatchArea;
	
    //Population structure
	vector<Patch> _Patch;

	//Resource consumption
    double  _FlucAmp;
    double  _FlucPeriod;
    double	_PatchInitialRes;
    double	_CoopEfficiency;
    double  _MaxIncreRate;//
	double	_MetabolicConsump;

	//Individual intake & fecundity
    double  _HalfConst;
    double  _ReproMax;
    double	_CostRate;

	//Survival rate
	double	_AgeStandard;
    double	_SurvRateUL;

	//Reproduction
	double	_Heritability;
    vector<double> _Offspring;
	
    //Random number generator
	RandGenerator RG;

    //For convenience of data computation
    double  Pi;
	bool	_End;
    int		_Span;
	int		_CutPoint;
    int     _LastNA;
	vector<double> _TSData;
    vector<double> _PatchReproData;
    vector<double> _IndivReproData;
	vector<double> _IndivDegData;
};



//Function that repeats the annual processes and that updates the state of population
void PopProcess::SystemUpdate(int Years){
	
	for ( int yy = 0 ; yy < Years ; yy++ ){
		
		//double FlucResource = _PatchInitialRes;
        double FlucResource = _PatchInitialRes + _FlucAmp*sin( 2*Pi*yy/_FlucPeriod );

        if ( yy == Years-1 ){
            _End = 1;
        }

        //Patch update
		for ( int i = 0 ; i < _PatchNum ; i++ ){
			_Patch[i].UpdateAge();
            _Patch[i].UpdateCoopNum();
            _Patch[i].UpdateCoopDegSum();
			_Patch[i].UpdateResource (FlucResource, _CoopEfficiency, _MaxIncreRate);
            _Patch[i].UpdateReproRate(_MetabolicConsump, _HalfConst, _ReproMax, _CostRate, _End);
        }
        
        if (_End){
            ComputePatchReproData();
        }
        
		OffspringDisperse();
		for ( int i = 0 ; i < _PatchNum ; i ++ ){
			_Patch[i].UpdateSurvRate(_SurvRateUL, _AgeStandard, _CostRate);
			_Patch[i].UpdateAlive();
		}
		
        //Compute and check population size
		_PopSize = 0;
		for ( int i = 0 ; i < _PatchNum ; i ++ ) _PopSize = _PopSize + _Patch[i].GetSize();
		if ( _PopSize == 0 ){
			ErrorHandling1();
            break;
		}
		if ( _PopSize > _PopSizeUL ){
			ErrorHandling2();
			break;
		}

        ComputePopCoopData();
        if ( _TotalCoopNum == 0 )   _LastNA = yy;

        //Compute time series data
		if ( yy%_Span == _Span-1 ){
			ComputeTimeSeriesData( yy/_Span, FlucResource );
		}
	}
}

class SetUpAndRun{
public:
	SetUpAndRun(double EnvironResDens, double bK, double CostRate, double CoopProp, int EndTime, int initial, int Round, int Span, double Amp, double Period, bool Social, bool Divided);
    void	ProcessShunt(int type);
	void	PopSummary();
	void	PopTimeSeries();
	void	FecunAndAltruistProp();
	void	PopFinalState();

private:
	int		_EndTime;
	int		_initial;
	int		_Round;
	int		_Span;
    int     _CutPoint;
    bool    _Social;
    bool    _Divided;
	double	_bK;
    double  _FlucAmp;
    double  _FlucPeriod;
	double	_CostRate;
	double	_CoopProp;
    double  _EnvironResDens;
};





int main(int arc, char *argv[]){

    int     type, EndTime, initial, Round, Span;
    bool    Social, Divided;
    double  EnvironResDens, bK, CostRate, CoopProp, Amp, Period;

    while ( cin >> type >> EnvironResDens >> bK >> CostRate >> CoopProp >> EndTime >> initial >> Round >> Span >> Amp >> Period >> Social >> Divided){
        SetUpAndRun SUR(EnvironResDens, bK, CostRate, CoopProp, EndTime, initial, Round, Span, Amp, Period, Social, Divided);
	    SUR.ProcessShunt(type);
    }
	
	return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

double StandardError(vector<double> data, double mean){
    double se = -1.0;
    if ( data.size() > 1 ){
        double temp = 0;
        for ( int i = 0; i < data.size(); i++ ){
            temp = temp + pow((data[i]-mean),2);
        }
        temp = temp/(data.size()-1);
	    se = sqrt(temp/data.size());
    } 
    return se;
}

RandGenerator::RandGenerator(){
	seed = chrono::system_clock::now().time_since_epoch().count();
	RAND_GENERATOR.seed(seed);
}

double RandGenerator::RUnif(){
	uniform_real_distribution<double> distribution(0.0,1.0);
	return distribution(RAND_GENERATOR);
}

double RandGenerator::Rand2(){
    uniform_int_distribution<int> distribution(1,10);
    double output = (double)distribution(RAND_GENERATOR)/10;
    return output;
}

bool RandGenerator::RBern(double p){
	bernoulli_distribution distribution(p);
	return distribution(RAND_GENERATOR);
}

int RandGenerator::SampleInt(int start, int end){
	uniform_int_distribution<int> distribution(start,end-1);
	return distribution(RAND_GENERATOR);
}

Individual::Individual(int PatchID, double CoopDeg, int Age){
	_PatchID        = PatchID;
	_Age            = Age;
	_CoopDeg        = CoopDeg;
	_ReproRate      = 0;

	if ( _CoopDeg == -1.0 )	_CoopDeg = RG.Rand2();
}

Patch::Patch(){
}

void Patch::AddIndividual(int PatchID, double CoopDeg, int Age){
	_Individual.push_back(Individual(PatchID, CoopDeg, Age));
}

void Patch::UpdateCoopDegSum(){
	_CoopDegSum = 0.0;
	for ( int i = 0 ; i < _Individual.size() ; i++ )
		_CoopDegSum = _CoopDegSum + _Individual[i].GetCoopDeg();
}

void Patch::UpdateCoopNum(){
	_CoopNum = 0;
    for ( int i = 0 ; i < _Individual.size() ; i++ ){
		if ( _Individual[i].GetCoopDeg() != 0.0 ){
			_CoopNum++;
    	}
	}
}

void Patch::UpdateAge(){
	for ( int i = 0 ; i < _Individual.size(); i++ )	_Individual[i].Aging();
}

void Patch::UpdateResource(double InitialResource, double CoopEfficiency, double MaxIncreRate){
    double tmp = _CoopDegSum * CoopEfficiency;
    _Resource = InitialResource*(1 + MaxIncreRate*tmp/(MaxIncreRate*InitialResource/2 + tmp));
}

void Patch::UpdateReproRate(double MetabolicConsump, double HalfConst, double ReproMax, double CostRate, bool End){
    //Compute individual intake
    double ResPerCap = _Resource/_Individual.size();
    double tmp1 = 0;
    if (ResPerCap > MetabolicConsump){
        tmp1 = (ResPerCap-MetabolicConsump) / (HalfConst+ ResPerCap-MetabolicConsump);
    }
	
	//Compute individual fecundity and total fecundity of patch //Definition of fecundity?
	_OffspringNum = 0;
	for ( int i = 0 ; i < _Individual.size() ; i++ ){
		double	tmp2 = ReproMax * tmp1 * ( 1.0 - CostRate*_Individual[i].GetCoopDeg() );
		double  tmp3 = floor(tmp2);
		double  tmpp = tmp2 - tmp3;
		int     IndivRepro;
		if (RG.RBern(tmpp)){
			IndivRepro = (int)(tmp3+1);
		}else{
			IndivRepro = (int)tmp3;
		}
		_OffspringNum = _OffspringNum + IndivRepro;
        if (End){
            _EveryIndivRepro.push_back(IndivRepro);
            _EveryIndivDeg.push_back(_Individual[i].GetCoopDeg());
        }
		_Individual[i].WriteReproRate(IndivRepro);
	}
}

void Patch::UpdateSurvRate(double SurvRateUL, double AgeStandard, double CostRate){
	for ( int i = 0 ; i < _Individual.size() ; i++ ){
		double IndivSurvRate = SurvRateUL*exp(-_Individual[i].GetAge()/AgeStandard );
		_Individual[i].WriteSurvRate(IndivSurvRate);
	}
}

void Patch::UpdateAlive(){
	for ( int i = 0 ; i < _Individual.size() ; i++ ){
		if ( !RG.RBern(_Individual[i].GetSurvRate()) ){
			_Individual.erase(_Individual.begin()+i);
			i--;
		}
	}
}

vector<double> Patch::GetOffspring(){
	vector<double> Offspring;
	for ( int i = 0 ; i < _Individual.size() ; i++ ){
		vector<double> tmp( _Individual[i].GetReproRate(), _Individual[i].GetCoopDeg() );
		Offspring.insert( Offspring.end(), tmp.begin(), tmp.end() );
	}
	return Offspring;
}



PopProcess::PopProcess(double CoopEfficiency, double EnvironResDens, double CostRate, double CoopProp, int Span, int CutPoint, double FlucAmp, double FlucPeriod, bool Social, bool Divided){
    //Population state
    _Divided              = Divided;
    _PopSizeUL            = 100000;
	_PopSize              = 300;

	if (Social)	    _CoopProp = CoopProp;
	else			_CoopProp = 0.0;
	if (Divided){
        _PatchNum   = 90;
        _PatchArea  = 1;
    }else{
        _PatchNum   = 1;
        _PatchArea  = 90;
    }            

    //Population structure
	_Patch = vector<Patch>(_PatchNum);

	//Resource consumption
    _FlucAmp              = FlucAmp;
    _FlucPeriod           = FlucPeriod;
	_PatchInitialRes      = EnvironResDens*_PatchArea;
	_CoopEfficiency       = CoopEfficiency;
    _MaxIncreRate         = 40;//
	_MetabolicConsump     = 1.0;

	//Individual intake & fecundity
    _HalfConst            = 2;
    _ReproMax             = 3.0;
    _CostRate          	  = CostRate;

	//Survival rate
	_AgeStandard          = 2.0;
    _SurvRateUL           = 0.7;

	//Reproduction
	if (Social)	_Heritability = 0.999;
	else		_Heritability = 1.0;

    //For convenience of data computation
    Pi              = 4*atan(1);
	_End            = 0;
    _Span           = Span;
	_CutPoint		= CutPoint;
    _LastNA         = -1;
	_TSData			= vector<double>(CutPoint *4, -1.0);
	_PatchReproData = vector<double>(_PatchNum*7, -1.0);
	
	//Initialize population
	for ( int i = 0 ; i < _PopSize ; i++ ){
		int PatchID = RG.SampleInt(0,_PatchNum);
		if ( i < _PopSize*_CoopProp ) _Patch[PatchID].AddIndividual(PatchID, -1.0, 0);
		else _Patch[PatchID].AddIndividual(PatchID, 0.0, 0);
	}
	ComputePopCoopData();
}

PopProcess::~PopProcess(){}

void PopProcess::OffspringDisperse(){
	_Offspring.clear();
	for ( int i = 0 ; i < _Patch.size() ; i++ ){
		vector<double> temp = _Patch[i].GetOffspring();
		_Offspring.insert( _Offspring.end(), temp.begin(), temp.end() );
	}
	for ( int i = 0 ; i < _Offspring.size() ; i++ ){
		int PatchID = RG.SampleInt(0,_PatchNum);
		if ( !RG.RBern(_Heritability) ){
			if ( _Offspring[i] != 0.0 ) _Offspring[i] = 0.0;
			else _Offspring[i] = RG.Rand2();
		}
        _Patch[PatchID].AddIndividual(PatchID, _Offspring[i], 0);
	}
}

void PopProcess::ErrorHandling1(){
	_TotalCoopNum   = 0;
	_CoopProp       = -1.0;
	_AverCoopDeg    = -1.0;
	_End            = 1;	
}

void PopProcess::ErrorHandling2(){
    _PopSize        = -1;
	_TotalCoopNum   = -1;
	_CoopProp       = -1.0;
	_AverCoopDeg    = -1.0;
	_End            = 1;
}

void PopProcess::ComputePopCoopData(){
	_TotalCoopNum = 0;
	_AverCoopDeg  = 0.0;
	for ( int i = 0 ; i < _Patch.size() ; i++ ){
		_Patch[i].UpdateCoopNum();
		_Patch[i].UpdateCoopDegSum();
		_TotalCoopNum = _TotalCoopNum + _Patch[i].GetCoopNum();
		_AverCoopDeg  = _AverCoopDeg  + _Patch[i].GetCoopDegSum();
	}
	_CoopProp    = double(_TotalCoopNum)/double(_PopSize);
    _AverCoopDeg = _AverCoopDeg/double(_PopSize);
}

void PopProcess::ComputePatchReproData(){
	//PatchSize\tcNum\tcProp\tcDeg\tcDegTotal\tProdTotal\tAverProd
    for ( int i = 0 ; i < _PatchNum ; i++ ){
		vector<double> tmp1 = _Patch[i].GetEveryIndivRepro();
		vector<double> tmp2 = _Patch[i].GetEveryIndivDeg();
		_IndivReproData.insert( _IndivReproData.end(), tmp1.begin(), tmp1.end() );
		_IndivDegData.insert( _IndivDegData.end(), tmp2.begin(), tmp2.end() );
        _Patch[i].UpdateCoopNum();
		_Patch[i].UpdateCoopDegSum();
        _PatchReproData[i*7]   = (double)_Patch[i].GetSize();
		_PatchReproData[i*7+1] = double(_Patch[i].GetCoopNum());
		_PatchReproData[i*7+4] = _Patch[i].GetCoopDegSum();
		_PatchReproData[i*7+5] = double(_Patch[i].GetOffspring().size());
        if ( _PatchReproData[i*7] == 0 ){
            _PatchReproData[i*7+2] = -1.0;
            _PatchReproData[i*7+3] = -1.0;
			_PatchReproData[i*7+6] = -1.0;
        }else{
            _PatchReproData[i*7+2] = double(_Patch[i].GetCoopNum())/double(_Patch[i].GetSize());
			_PatchReproData[i*7+3] = _Patch[i].GetCoopDegSum()/double(_Patch[i].GetSize());
			_PatchReproData[i*7+6] = double(_Patch[i].GetOffspring().size())/double(_Patch[i].GetSize());
        }
    }
}

void PopProcess::ComputeTimeSeriesData(int index, double FlucRes){
	_TSData[index]				= _PopSize;
	_TSData[index+1*_CutPoint]	= _CoopProp;
	_TSData[index+2*_CutPoint]	= _AverCoopDeg;
	_TSData[index+3*_CutPoint]  = FlucRes;
}

SetUpAndRun::SetUpAndRun(double EnvironResDens, double bK, double CostRate, double CoopProp, int EndTime, int initial, int Round, int Span, double Amp, double Period, bool Social, bool Divided){
	_EnvironResDens 	= EnvironResDens;
    _FlucAmp            = Amp;
    _FlucPeriod         = Period;
	_EndTime			= EndTime;
	_initial			= initial;
	_Round				= Round;
	_Span				= Span;
    _CutPoint           = EndTime/Span;
	_bK					= bK;
	_CostRate			= CostRate;
	_Social			    = Social;
    _Divided            = Divided;
	_CoopProp			= CoopProp;
}

void SetUpAndRun::ProcessShunt(int type){
	switch (type)
	{
		case 1:{
			PopSummary();
		}break;
		case 2:{
			PopTimeSeries();
		}break;
		case 3:{
			FecunAndAltruistProp();
		}break;
		case 4:{
			PopFinalState();
		}break;
	}
}

void SetUpAndRun::PopSummary(){
	//cout << "ResourceLevel\tCooperationEfficiency\tCostRate\tTotalRound\tPopSize\tpSizeSe\tAltruistProp\taPropSe\tAverAltruismDeg\taDegSe\n"; 
    vector<int>                 count( 3     , 0  );
	vector<double>              total( 3     , 0.0);
	vector<double>              row  ( _Round, 0.0);
	vector< vector<double> >    table( 3     , row);
	for ( int i = 0; i < _Round; i++ ){
		PopProcess PP(_bK, _EnvironResDens, _CostRate, _CoopProp, _Span, _CutPoint, _FlucAmp, _FlucPeriod, _Social, _Divided);
        PP.SystemUpdate(_EndTime);
		double final_state[3] = {(double)PP.GetPopSize(), PP.GetCoopProp(), PP.GetAverCoopDeg()};
		for ( int j = 0; j < 3; j++ ){
			if ( final_state[j] != -1.0 ){
				table[j][count[j]] = final_state[j];
				total[j] = total[j] + final_state[j];
				count[j]++;
			}
		}
	}
	cout << _EnvironResDens << '\t' << _bK << '\t' << _CostRate << '\t' << _Round;
	for ( int j = 0; j < 3; j++ ){
		if ( count[j] != 0 ){
			double mean = total[j]/count[j];
			vector<double> tmp;
			vector<double>::iterator it;
			it = table[j].begin();
			tmp.assign(it, it+count[j]);
			cout << '\t' << mean << '\t' << StandardError(tmp, mean);
		}else{
			cout << "\t-1\t-1";
		}
	}
	cout << '\n';
}

void SetUpAndRun::PopTimeSeries(){

	//cout << "ResourceLevel\tCooperationEfficiency\tCostRate\tRound\tTime\tPopSize\tAltruistProp\tAverAltruismDeg\tIC\n";
	
	for ( int i = 0; i < _Round; i++ ){
		PopProcess PP(_bK, _EnvironResDens, _CostRate, _CoopProp, _Span, _CutPoint, _FlucAmp, _FlucPeriod, _Social, _Divided);
        cout << _EnvironResDens << '\t' << _bK << '\t' << _CostRate << '\t' << _initial+i << "\t0\t" << _FlucAmp << '\t' << _FlucPeriod << '\t' \
             << PP.GetPopSize() << '\t' << PP.GetCoopProp() << '\t' << PP.GetAverCoopDeg() << '\t' << _CoopProp << '\n';

		PP.SystemUpdate(_EndTime);
		vector<double> TSData_tmp = PP.GetTS();
		for ( int j = 0; j < _CutPoint; j++ ){
            cout << TSData_tmp[j+_CutPoint*3] << '\t' << _bK << '\t' << _CostRate << '\t' << _initial+i << '\t' << (j+1)*_Span << '\t' << _FlucAmp << '\t' << _FlucPeriod;
			for ( int k = 0; k < 3; k++ )   cout << '\t' << TSData_tmp[j+_CutPoint*k];
            cout << '\t' << _CoopProp << '\n';
		}
	}
}

void SetUpAndRun::FecunAndAltruistProp(){
	//cout << "ResDensity\tCooperationEfficiency\tCostRate\tEndTime\tRound\tPatchID\tIndivDeg\tIndivProd\tPatchSize\tcNum\tcProp\tcDeg\tcDegTotal\tProdTotal\tAverProd\n";
	
	for ( int i = 0; i < _Round; i++ ){
        PopProcess PP(_bK, _EnvironResDens, _CostRate, _CoopProp, _Span, _CutPoint, _FlucAmp, _FlucPeriod, _Social, _Divided);
	    PP.SystemUpdate(_EndTime);
        int base = 0;
	    vector<double> PatchData_tmp  = PP.GetPatchReproData();
        vector<double> IndivFecun_tmp = PP.GetIndivReproData();
		vector<double> IndivDeg_tmp   = PP.GetIndivDegData();
        for ( int j = 0; j < 90; j++ ){
            for ( int k = 0; k < PatchData_tmp[j*7]; k++ ){
                cout << _EnvironResDens << '\t' << _bK << '\t' << _CostRate << '\t' << _EndTime << '\t' \
					 << _initial+i << '\t' << j << '\t' << IndivDeg_tmp[base+k] << '\t' << IndivFecun_tmp[base+k] \
					 << '\t' << PatchData_tmp[j*7] << '\t' << PatchData_tmp[j*7+1] << '\t' \
					 << PatchData_tmp[j*7+2] << '\t' << PatchData_tmp[j*7+3] << '\t' << PatchData_tmp[j*7+4] \
					 << '\t' << PatchData_tmp[j*7+5] << '\t' << PatchData_tmp[j*7+6] << '\n';
            }
            base = base + PatchData_tmp[j*7];
	    }
    }
}

void SetUpAndRun::PopFinalState(){
	//cout << "ResourceLevel\tCooperationEfficiency\tCostRate\tState\tRound\tPopSize\tAltruistProp\tAverAltruismDeg\tLastNA\tIC\n";
	for ( int i = 0; i < _Round; i++ ){
		PopProcess PP(_bK, _EnvironResDens, _CostRate, _CoopProp, _Span, _CutPoint, _FlucAmp, _FlucPeriod, _Social, _Divided);
		PP.SystemUpdate(_EndTime);
		cout << _EnvironResDens << '\t' << _bK << '\t' << _CostRate;
        if ( _Social )	cout << "\tSocial\t";
		else			cout << "\tNonsocial\t";
		cout << _initial+i << '\t' << PP.GetPopSize() << '\t' << PP.GetCoopProp() << '\t' \
			 << PP.GetAverCoopDeg() << '\t' << PP.GetLastNA() << '\t' << _CoopProp << '\n';
	}
}
