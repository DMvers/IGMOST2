#include <regex>
#include <string>
#include "BactParam.h"

/* These functions assist the class BactParam (AGORAparam) to read in the SBML files.
 */

bool hasPrefix(std::string const &fullString, std::string const &prefix)
{
    // Check if a string starts with a certain substring.
    if (fullString.length() >= prefix.length())
    {
        return (0 == fullString.compare(0, prefix.length(), prefix));
    }
    else
    {
        return false;
    }
}

bool hasEnding(std::string const &fullString, std::string const &ending)
{
    // Check if a string ends with a certain substring.
    if (fullString.length() >= ending.length())
    {
        return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    }
    else
    {
        return false;
    }
}

std::string make_universal_id(std::string id)
{
    // Take an id string from a SBML model and make it into an universal id that is the same for all models.
    // These Regex functions are famously slow, but we only run them during initialisation, so it's not a problem for us
    std::regex reg1("(__|_DASH_)");
    std::regex reg2("^M_");
    std::regex reg3("_[epci]$");
    std::regex reg4("_LSQBKT[ec]_RSQBKT$");
    std::regex reg5("^bigg_");
    std::regex reg6("_bm$");
    std::regex reg7("_91_[ec]_93_$");
    std::regex reg8("_40_[ec]_41_$");
    id = std::regex_replace(id, reg1, "_");
    id = std::regex_replace(id, reg2, "");
    id = std::regex_replace(id, reg3, "");
    id = std::regex_replace(id, reg4, "");
    id = std::regex_replace(id, reg5, "");
    id = std::regex_replace(id, reg6, "");
    id = std::regex_replace(id, reg7, "");
    id = std::regex_replace(id, reg8, "");
    return id;
}

int getCarbon(std::string notes)
{
    // Get amount of carbon atoms in a metabolite.
    std::string formula;
    std::string x = "fbc:chemicalFormula=\"";
    int index;
    index = notes.find(x);
    if (index == -1)
        return 0;
    formula = notes.substr(index + x.length());
    index = formula.find("\"");
    formula = formula.substr(0, index);
    std::regex reg("(.*C)([0-9]*)([A-Z].*)");
    std::string number = std::regex_replace(formula, reg, "$2");
    if (number.compare(formula) == 0)
        return 0;
    else if (number.compare("") == 0)
        return 1;
    return atoi(number.c_str());
}

// Used to get the full English name of a metabolite
std::string getFullName(std::string notes)
{
    std::string formula;
    std::string x = "name=\"";
    int index;
    index = notes.find(x);
    if (index == -1)
        return "?";
    formula = notes.substr(index + x.length());
    index = formula.find("\"");
    formula = formula.substr(0, index);
    return formula;
}

std::string getKeggName(std::string notes)
{
    std::string formula;
    std::string x = "kegg.compound/"; // for finding the name
    // std::string x = "rdf:resource=\"";
    int index;
    index = notes.find(x);
    if (index == -1)
        return "Unknown";
    formula = notes.substr(index + x.length());
    index = formula.find("\""); // for finding the end of the name
    // index = formula.find("<"); //end of chebiid
    formula = formula.substr(0, index);
    // std::cout<<formula;
    return formula;
    // std::cout<<"formula is " <<formula <<std::endl;
}

std::string getChebiName(std::string notes)
{
    std::string formula;
    std::string x = "ChEBIID:";
    int index;
    index = notes.find(x);
    if (index == -1)
        return "Unknown";
    formula = notes.substr(index + x.length());
    index = formula.find("<"); // end of chebiid
    formula = formula.substr(0, index);
    return formula;
}

// Used to get the number of hydrogen atoms in a molecule (metabolite)
int getHydrogen(std::string notes)
{
    std::string formula;
    std::string x = "fbc:chemicalFormula=\"";
    int index;
    index = notes.find(x);
    if (index == -1)
        return 0;
    formula = notes.substr(index + x.length());
    index = formula.find("\"");
    formula = formula.substr(0, index);
    // Get hydrogen
    std::regex getH("H[[:digit:]]*");
    std::smatch match;
    std::string result = "nomatch";
    if (std::regex_search(formula, match, getH))
    {
        result = match.str();
        int length = result.length();
        if (length == 1)
        {
            return 1;
        }
        else
        {
            return stoi(result.substr(1, length));
        }
    }
    return 0;
}

// Used to get the number of nitrogen atoms in a molecule (metabolite)
int getNitrogen(std::string notes)
{
    std::string formula;
    std::string x = "fbc:chemicalFormula=\"";
    int index;
    index = notes.find(x);
    if (index == -1)
        return 0;
    formula = notes.substr(index + x.length());
    index = formula.find("\"");
    formula = formula.substr(0, index);
    // Get nitrogen
    std::regex getN("N[[:digit:]]*");
    std::smatch match;
    std::string result = "nomatch";
    if (std::regex_search(formula, match, getN))
    {
        result = match.str();
        int length = result.length();
        if (length == 1)
        {
            return 1;
        }
        else
        {
            return stoi(result.substr(1, length));
        }
    }

    return 0;
}

// Used to get the number of nitrogen atoms in a molecule (metabolite)
int getOxygen(std::string notes)
{
    std::string formula;
    std::string x = "fbc:chemicalFormula=\"";
    int index;
    index = notes.find(x);
    if (index == -1)
        return 0;
    formula = notes.substr(index + x.length());
    index = formula.find("\"");
    formula = formula.substr(0, index);
    // Get nitrogen
    std::regex getO("O[[:digit:]]*");
    std::smatch match;
    std::string result = "nomatch";
    if (std::regex_search(formula, match, getO))
    {
        result = match.str();
        int length = result.length();
        if (length == 1)
        {
            return 1;
        }
        else
        {
            return stoi(result.substr(1, length));
        }
    }

    return 0;
}

int contains_pair(QVector<int> v1, QVector<int> v2, int x1, int x2)
{
    // Find the index of a pair in two vectors.
    for (int a = 0; a < v1.size(); ++a)
        if (v1[a] == x1 && v2[a] == x2)
            return a;
    return -1;
}
