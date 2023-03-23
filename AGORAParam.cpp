#include <iostream>
#include "BactParam.h"
#include <regex>
#include <sstream>
#include "parameter.h"

using namespace std;
extern Par par;
extern std::string simulationdir;

AGORAParam::AGORAParam(std::string file, int i, double fw, std::map<std::string, double> gibbsmap)
{
    /* Constructs an object on the basis of an SBML file (file) and a value for the flux weight (fw).
     * Reads in the document and passes it to createDicts and addBacterium.
     */
    indexnr = i;
    FW = fw;
    SBMLReader reader;
    Model *model;
    SBMLDocument *doc = reader.readSBML(file);
    if (doc->getNumErrors() > 0)
    {
        std::cout << doc->getError(0)->getMessage() << std::endl
                  << "And " << (doc->getNumErrors()) - 1 << " more errors in the following model: " << file << std::endl;
    }
    model = doc->getModel();
    createDicts(model, gibbsmap);
    addBacterium(model);
    // std::cout<<gibbsmap["gal"];
}

AGORAParam::~AGORAParam()
{
    /* Destructor, does nothing in this case. */
}

void AGORAParam::createDicts(Model *m, std::map<std::string, double> gibbsmap)
{
    /* This method fills the maps/dictionaries react_dict, metab_dict and exchange_dict
     * and adds the indexes of the exchange reactions to the vector exchange_list.
     */
    numreactions = 0;
    nummetabolites = 0;
    int nmetabolites = m->getNumSpecies();
    ListOfSpecies *lmetabolites = m->getListOfSpecies();
    bool printexchange = false;
    // Enable the printing, to console and file, of all exchanged metabolites, if desired
    if (printexchange)
    {
        std::string metabnamefilename;
        std::string metabidfilename;
        std::string metabchebifilename;
        std::string metabkeggfilename;
        metabnamefilename = simulationdir + "/exchangenamelist" + m->getName() + ".csv";
        metabidfilename = simulationdir + "/exchangeidlist" + m->getName() + ".csv";
        metabchebifilename = simulationdir + "/exchangechebilist" + m->getName() + ".csv";
        metabkeggfilename = simulationdir + "/exchangekegglist" + m->getName() + ".csv";
        EXCHANGENAMEFILE.open(metabnamefilename);
        EXCHANGEIDFILE.open(metabidfilename);
        EXCHANGEKEGGFILE.open(metabkeggfilename);
        EXCHANGECHEBIFILE.open(metabchebifilename);
        std::cout << metabnamefilename;
    }

    // make the metabolite dictionary - each is assigned a number
    for (int j = 0; j < nmetabolites; ++j)
    {
        std::string id = lmetabolites->get(j)->getId();
        if (not metab_dict.count(id))
        {
            std::string compartment = lmetabolites->get(j)->getCompartment();
            metab_dict[id] = ++nummetabolites;
        }
    }
    int nreactions = m->getNumReactions();
    ListOfReactions *lreactions = m->getListOfReactions();
    for (int j = 0; j < nreactions; ++j)
    {
        std::string id = lreactions->get(j)->getId();
        if (not react_dict.count(id))
        {
            react_dict[id] = ++numreactions;
            if (hasPrefix(id, "R_EX") || hasPrefix(id, "R_DM_") || hasPrefix(id, "R_sink_")) // These are the three types of reaction that can add or remove metabolites from the FBA system
            {
                std::string univ_id = make_universal_id(lreactions->get(j)->getReactant(0)->getSpecies());
                if (printexchange)
                {
                    XMLNode *speciesnode = m->getSpecies(univ_id)->toXMLNode(); //->toXMLString();
                    std::string notes = speciesnode->toXMLString();
                    EXCHANGENAMEFILE << getFullName(notes) << std::endl;
                    EXCHANGEIDFILE << univ_id << std::endl;
                    EXCHANGEKEGGFILE << getKeggName(notes) << std::endl;
                    EXCHANGECHEBIFILE << getChebiName(notes) << std::endl;
                    delete speciesnode;
                }
                // Get gibbs energy
                if (gibbsmap.find(univ_id) != gibbsmap.end())
                {
                    exchange_dict[univ_id].energy = gibbsmap[univ_id];
                }
                else
                {
                    exchange_dict[univ_id].energy = 0.0;
                }
                XMLNode *speciesnode = m->getSpecies(lreactions->get(j)->getReactant(0)->getSpecies())->toXMLNode(); //->toXMLString();
                std::string notes = speciesnode->toXMLString();
                exchange_dict[univ_id].carbon = getCarbon(notes);
                exchange_dict[univ_id].hydrogen = getHydrogen(notes);
                exchange_dict[univ_id].nitrogen = getNitrogen(notes);
                exchange_dict[univ_id].oxygen = getOxygen(notes);
                exchange_dict[univ_id].reaction_out = numreactions;
                exchange_dict[univ_id].reaction_in = numreactions + 1;
                delete speciesnode;
            }
        }
        if (lreactions->get(j)->getReversible() || hasPrefix(id, "R_EX") || hasPrefix(id, "R_DM_") || hasPrefix(id, "R_sink_"))
        {
            std::string r_id = id + "_r";
            if (not react_dict.count(r_id))
            {
                react_dict[r_id] = ++numreactions;
            }
        }
    }
    for (auto const &x : exchange_dict)
    {
        exchange_list.append(x.second.reaction_in);
        exchange_list.append(x.second.reaction_out);
    }
    // Close the files we printed to, if we opened them
    if (printexchange)
    {
        EXCHANGENAMEFILE.flush();
        EXCHANGEIDFILE.flush();
        EXCHANGEKEGGFILE.flush();
        EXCHANGECHEBIFILE.flush();
    }
}

void AGORAParam::addReaction(Reaction *react)
{
    /*  This method adds the reactions of react_dict to the vectors genome, lbs, ubs, fluxs and objs */
    int j = react_dict.at(react->getId());
    genome[j - 1] = 1;
    lbs[j - 1] = 0;
    ubs[j - 1] = 10000000;
    fluxs[j - 1] = 0;
    objs[j - 1] = 0;
    if (react->getReversible())
    {
        j = react_dict.at(react->getId() + "_r");
        genome[j - 1] = 1;
        lbs[j - 1] = 0;
        ubs[j - 1] = 10000000;
        fluxs[j - 1] = 0;
        objs[j - 1] = 0;
    }
}

std::vector<double> AGORAParam::addReactants(Model *m, Reaction *react)
{
    /* This method adds the stoichiometric values of the reactants (<0) to the matrix S. */
    int nreactants = react->getNumReactants();
    int j1 = react_dict.at(react->getId());
    std::vector<double> elements = {0, 0, 0, 0};
    for (int y = 0; y < nreactants; ++y)
    {
        std::string id = react->getReactant(y)->getSpecies();
        if (not metab_dict.count(id))
            continue;
        int i = metab_dict.at(id);
        double stoich = react->getReactant(y)->getStoichiometry();
        int a = contains_pair(iRow, jCol, i, j1);
        XMLNode *speciesnode = m->getSpecies(id)->toXMLNode(); //->toXMLString();
        std::string notes = speciesnode->toXMLString();
        elements.at(0) += (getCarbon(notes) * stoich);
        elements.at(1) += (getHydrogen(notes) * stoich);
        elements.at(2) += (getNitrogen(notes) * stoich);
        elements.at(3) += (getOxygen(notes) * stoich);
        delete speciesnode;
        // std::cout<<getCarbon(m,id)<<" "<<stoich<<" "<<id<<" "<<react->getId()<<std::endl;
        if (-1 != a)
        {
            SparseS[a] += stoich * -1.0;
        }
        else
        {
            iRow.append(i);
            jCol.append(j1);
            SparseS.append(stoich * -1.0);
        }

        if (react->getReversible())
        {
            int j2 = react_dict.at(react->getId() + "_r");
            a = contains_pair(iRow, jCol, i, j2);
            if (-1 != a)
            {
                SparseS[a] += stoich;
            }
            else
            {
                iRow.append(i);
                jCol.append(j2);
                SparseS.append(stoich);
            }
        }
    }
    // std::cout<<elements[0] << std::endl;
    return (elements);
}

std::vector<double> AGORAParam::addProducts(Model *m, Reaction *react)
{
    /* This method adds the stoichiometric values of the products (>0) to the matrix S. */
    int nproducts = react->getNumProducts();
    int j1 = react_dict.at(react->getId());
    std::vector<double> elements = {0, 0, 0, 0};
    for (int y = 0; y < nproducts; ++y)
    {
        std::string id = react->getProduct(y)->getSpecies();
        if (not metab_dict.count(id))
            continue;
        int i = metab_dict.at(id);
        double stoich = react->getProduct(y)->getStoichiometry();

        int a = contains_pair(iRow, jCol, i, j1);
        XMLNode *speciesnode = m->getSpecies(id)->toXMLNode(); //->toXMLString();
        std::string notes = speciesnode->toXMLString();
        elements.at(0) += (getCarbon(notes) * stoich);
        elements.at(1) += (getHydrogen(notes) * stoich);
        elements.at(2) += (getNitrogen(notes) * stoich);
        elements.at(3) += (getOxygen(notes) * stoich);
        delete speciesnode;
        if (-1 != a)
        {
            SparseS[a] += stoich;
        }
        else
        {
            iRow.append(i);
            jCol.append(j1);
            SparseS.append(stoich);
        }
        if (react->getReversible())
        {
            int j2 = react_dict.at(react->getId() + "_r");
            a = contains_pair(iRow, jCol, i, j2);
            if (-1 != a)
            {
                SparseS[a] += stoich * -1.0;
            }
            else
            {
                iRow.append(i);
                jCol.append(j2);
                SparseS.append(stoich * -1.0);
            }
        }
    }
    // std::cout<<elements[0] << std::endl;
    return (elements);
}

void AGORAParam::setBiomassReaction(Model *m)
{
    /* This method creates the standard biomass reaction used to replace the biomass reaction present in the SBML model. */
    // alternatively, it finds an annotates the default biomass reaction to work nicely with our checks on N,C,H balance

    // We search through all keys for one that starts with biomass
    // Inefficient, but only used during startup, so not unaceptably so
    int j;
    std::string biomassnamestringed;
    for (int i = 0; i < 1000; i++)
    {
        std::ostringstream biomassname;
        biomassname << "R_biomass";
        if (i < 10)
        {
            biomassname << 0;
        }
        if (i < 100)
        {
            biomassname << 0;
        }
        biomassname << i;

        biomassnamestringed = biomassname.str();
        if (react_dict.count(biomassnamestringed) > 0)
        {
            j = react_dict.at(biomassnamestringed);
            break;
        }
    }
    bm_reaction = j;
    objs[j - 1] = 1.0; // Set the function we found as the objective function

    if (par.knockoutxfp)
    {
        std::set<string> bifids = {"Bbre", "Binf", "Blong"};
        if (bifids.find(m->getName()) != bifids.end())
        {
            std::cout << "Removing reactions from " << m->getName() << std::endl;
            remove_reaction("R_PKL");
            remove_reaction("R_F6PE4PL");
        }
    }

    if (par.knockoutbiflactate || par.knockoutbifandecollactate)
    {
        std::set<string> bifids = {"Bbre", "Binf", "Blong"};
        if (bifids.find(m->getName()) != bifids.end())
        {
            std::cout << "Removing reactions from " << m->getName() << std::endl;
            remove_reaction("R_LDH_L");
        }
    }

    if (par.knockoutbiflactateup)
    {
        std::set<string> bifids = {"Bbre", "Binf", "Blong"};
        if (bifids.find(m->getName()) != bifids.end())
        {
            std::cout << "Removing reactions from " << m->getName() << std::endl;
            remove_reaction("R_LDH_L_r");
        }
    }

    if (par.knockoutecollactate || par.knockoutbifandecollactate)
    {
        std::set<string> bifids = {"Ecol"};
        if (bifids.find(m->getName()) != bifids.end())
        {
            std::cout << "Removing reactions from " << m->getName() << std::endl;
            remove_reaction("R_LDH_L");
        }
    }

     if(par.knockoutbutyro12ppd){
        std::set<string> butyrogenics = {"Cbut","Rinu","Ehal"};
            if(butyrogenics.find(m->getName())!=butyrogenics.end()){
            std::cout<<"Removing reactions from " <<m->getName()<<std::endl;
            remove_reaction("R_12PPDt");
        }
     }
    if(par.knockoutbutyrolactate){
        std::set<string> butyrogenics = {"Cbut","Rinu","Ehal"};
            if(butyrogenics.find(m->getName())!=butyrogenics.end()){
            std::cout<<"Removing reactions from " <<m->getName()<<std::endl;
            remove_reaction("R_EX_lac_L__40__e__41___r"); //TODO
            remove_reaction("R_EX_lac_D__40__e__41___r"); //TODO
        }
    }
    if(par.knockoutbutyrolcts){
        std::set<string> butyrogenics = {"Cbut","Rinu","Ehal"};
            if(butyrogenics.find(m->getName())!=butyrogenics.end()){
            std::cout<<"Removing reactions from " <<m->getName()<<std::endl;
            remove_reaction("R_EX_lcts__40__e__41___r"); //TODO
        }
    
    }
    

    if (par.biomassreplace) // Replace the old function in place with our own, simpler function
    {
        // Remove old reaction
        while (true)
        {
            int x = jCol.indexOf(j);
            if (x == -1)
                break;
            iRow.remove(x);
            jCol.remove(x);
            SparseS.remove(x);
        }

        // New Biomass reaction: ATP + H2O -> ADP +pi + h + (weightless)"biomass"
        iRow.append(metab_dict.at("M_h2o__91__c__93__"));
        jCol.append(j);
        SparseS.append(-50);
        iRow.append(metab_dict.at("M_atp__91__c__93__"));
        jCol.append(j);
        SparseS.append(-50);
        iRow.append(metab_dict.at("M_adp__91__c__93__"));
        jCol.append(j);
        SparseS.append(50);
        iRow.append(metab_dict.at("M_pi__91__c__93__"));
        jCol.append(j);
        SparseS.append(50);
        iRow.append(metab_dict.at("M_h__91__c__93__"));
        jCol.append(j);
        SparseS.append(50); // This should be 50
    }
    else // If not replacing biomass, we should annotate the old biomass version with the right value for Carbon, Hydrogen en Nitrogen
    {
        double carbon = 0;
        double nitrogen = 0;
        double hydrogen = 0;
        double oxygen = 0;
        int x = jCol.indexOf(j);
        while (x <= jCol.lastIndexOf(j))
        {
            // if (x == -1) break;
            std::string metabid;
            for (auto it = metab_dict.begin(); it != metab_dict.end(); ++it)
            {

                if (it->second == iRow[x])
                {
                    // Found metabolite, now to find the strength
                    metabid = it->first;
                    break;
                }
            }
            double strength = SparseS[x];                               // Finding strength is easy, can just be done by index number
            XMLNode *speciesnode = m->getSpecies(metabid)->toXMLNode(); //->toXMLString();
            std::string notes = speciesnode->toXMLString();
            carbon -= strength * getCarbon(notes);
            nitrogen -= strength * getNitrogen(notes);
            hydrogen -= strength * getHydrogen(notes);
            oxygen -= strength * getOxygen(notes);
            delete speciesnode;
            x++;
        }
        // Save the biomass weight using the species name. Easier than using the biomass name, because we work with species names in a lot of places already
        exchange_dict[m->getName()].carbon = carbon;
        exchange_dict[m->getName()].nitrogen = nitrogen;
        exchange_dict[m->getName()].hydrogen = hydrogen;
        exchange_dict[m->getName()].oxygen = oxygen;
        std::cout << "This model biomass has the following C, N, H,: " << carbon << " " << nitrogen << " " << hydrogen << std::endl;
    }
}

void AGORAParam::remove_reaction(std::string reaction)
{
    // This can be used to remove ethanol export
    if (react_dict.count(reaction) > 0)
    {
        std::cout << "Removing " << reaction << std::endl;
        int j = react_dict.at(reaction);
        // Remove reaction
        while (true)
        {
            int x = jCol.indexOf(j);
            if (x == -1)
                break;
            iRow.remove(x);
            jCol.remove(x);
            SparseS.remove(x);
        }
    }
    else
    {
        std::cout << "Reaction " << reaction << " not found" << std::endl;
    }
}

void AGORAParam::add_FWs()
{
    /* This method adds the flux weights (one value for all reactions in this case) to the matrix S. */
    jCol.resize(numstoich + numreactions);
    iRow.resize(numstoich + numreactions);
    SparseS.resize(numstoich + numreactions);
    for (int i = numstoich; i < numstoich + numreactions; i++)
    {
        std::string id = "blank";
        for (auto it = react_dict.begin(); it != react_dict.end(); ++it)
            if (it->second == (i - numstoich))
            {
                id = it->first;
                // std::cout<<id<<std::endl;
            }
        iRow[i] = nummetabolites + 1;
        jCol[i] = i - numstoich + 1;
        if (hasPrefix(id, "R_EX") || hasPrefix(id, "R_DM_") || hasPrefix(id, "R_sink_")) // These are the three types of reaction that can add or remove metabolites from the FBA system
        {
            // std::cout<<id<<std::endl;
            SparseS[i] = FW; // FW;
        }
        else
        {
            SparseS[i] = FW;
        }
        // lreactions->get(i)->getReactant(0)->getSpecies()
    }
}

void AGORAParam::addBacterium(Model *m)
{
    /* This method calls all the methods to load one bacterium species into the class. */
    speciesName = m->getName();
    genome.resize(numreactions);
    used.resize(numreactions);
    lbs.resize(numreactions);
    ubs.resize(numreactions);
    fluxs.resize(numreactions);
    objs.resize(numreactions);

    iRow.append(0);
    jCol.append(0);
    SparseS.append(0);

    int nreactions = m->getNumReactions();
    ListOfReactions *lreactions = m->getListOfReactions();
    Reaction *reaction;
    for (int x = 0; x < nreactions; ++x)
    {
        reaction = lreactions->get(x);
        if (react_dict.count(reaction->getId()))
        {
            addReaction(reaction);
            std::vector<double> inputs;
            std::vector<double> outputs;
            inputs = addReactants(m, reaction);
            outputs = addProducts(m, reaction);
            if (hasPrefix(reaction->getId(), "R_EX") || hasPrefix(reaction->getId(), "R_DM_") || hasPrefix(reaction->getId(), "R_biomass") || hasPrefix(reaction->getId(), "R_sink_"))
            {
                continue;
            }
            for (int i = 0; i < inputs.size(); i++)
            {
                if (inputs[i] != outputs[i])
                {
                    std::cout << "Error in reaction"
                              << " " << speciesName << " " << reaction->getId() << " " << inputs[i] << " " << outputs[i] << " " << i << std::endl;
                }
            }
        }
    }
    setBiomassReaction(m);

    numstoich = SparseS.size();
    add_FWs();

    lp = glp_create_prob();
    params = (glp_smcp *)malloc(sizeof(glp_smcp));
    glp_init_smcp(params);

    glp_add_rows(lp, nummetabolites + 1);
    glp_add_cols(lp, numreactions);

    // As metabolites are fixed, only do this at initialisation
    for (int i = 1; i <= nummetabolites + 1; i++)
    {
        if (i <= nummetabolites)
            glp_set_row_bnds(lp, i, GLP_FX, 0, 0);
        else if (i == nummetabolites + 1)
            glp_set_row_bnds(lp, i, GLP_DB, 0, 0.2); // This 0.2 is used in reaching a flux limit
    }

    // Load the matrix into GLPK:
    glp_load_matrix(lp, numstoich + numreactions - 1, iRow.data(), jCol.data(), SparseS.data());

    params->presolve = false; // whether presolve is on. presolving discards the previous found solution
    params->msg_lev = 0;
    params->tm_lim = 200000;
}
