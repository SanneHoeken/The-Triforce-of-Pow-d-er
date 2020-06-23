# Triforce of Pow(d)er

Dit project valt onder het vak Programmeertheorie dat onderdeel is van de minor Programmeren aan de Universiteit van Amsterdam.
Juni 2020.

### Case: Protein Pow(d)er

Eiwitten zijn lange strengen van aminozuren die aan de basis staan van veel belangrijke processen in het menselijk lichaam. Bij het maken van eiwitten wordt, na het synthetiseren van een keten van aminzuren, een eiwit in een 3D-structuur gevouwen. Een specifieke vouwing is bepalendvoor de functie van een eiwit. Het is daarom van groot belang dat een eiwit zich goed opvouwt; verkeerd gevouwen eiwitten kunnen namelijk resulteren in kanker, Alzheimer of taaislijmziekte. Een belangrijke factor in de vouwing van eiwitten is de stabiliteit van het eiwit. Aantrekkingskrachten tussen hydrofobe aminozuren (H) en Cysteine aminozuren(C), in tegenstelling tot polaire aminozuren (P), zorgen ervoor dat tegenover elkaar liggende H of C'tjes een bond vormen die het eiwit stabieler maakt. Het doel van deze case is een gegeven eiwit zo stabiel mogelijk op te vouwen (dus zo veel mogelijk H-bonds te genereren). 

In deze case gebruiken we hiervoor een model waarin de aminozuren in een eiwit op een twee- of driedimensionaal grid geplaatst worden. Elk aminozuur komt op een gridpunt te liggen en een volgend aminozuur ligt op één van de aangrenzende gridpunten met hoeken van 90 graden. Als twee H’s naast elkaar op het grid liggen krijgt het totale eiwit een -1 op de score. Als twee Cysteine-aminozuren naast elkaar liggen krijgt het eiwit -5 op de score. Tussen C’s en H’s is de score -1, en met P’s is er geen bindingseffect, dus score nul. Hoe lager de score, hoe stabieler het eiwit.

## Aan de slag

### Vereisten

Deze codebase is volledig geschreven in Python 3.7. In requirements.txt staan alle benodigde packages om de code succesvol te draaien. Deze zijn gemakkelijk te installeren via pip dmv. de volgende instructie:

```
pip install -r requirements.txt
```

Of via conda:

```
conda install --file requirements.txt
```

### Gebruik

Een voorbeeldje kan gerund worden door aanroepen van:

```
python main.py
```

Het gebruiksvriendelijke programma zal je vragen een eiwit-string in te voeren (in hoofdletters en zonder spaties), en vervolgens de gewenste dimensionaliteit, het gewenste algoritme en waar nodig de gewenste parameters in te voeren. Na de uitvoering van het algoritme word je gevraagd een path naar een csv-file in te voeren waarin de output zal worden opgeslagen. Ten slotte krijg je een plot van het resultaat te zien.

### Structuur

De hierop volgende lijst beschrijft de belangrijkste mappen en files in het project, en waar je ze kan vinden:

- **/code**: bevat alle code van dit project
  - **/code/algorithms**: bevat de code voor algoritmes
    - **/code/algorithms/help_methods**: bevat overige functies
  - **/code/classes**: bevat de benodigde classes voor deze case
  - **/code/visualisation**: bevat de code voor de visualisatie
- **/data**: bevat een tekstbestand met eiwit-strings die als input zijn te gebruiken

## Algoritmes

1. **Random.**
  Het Random algoritme vouwt voor elke iteratie een eiwit door willekeurige mogelijke vouwwaarden één voor één toe te wijzen aan elk aminozuur. Als een aminozuur eindigt in een doodlopende ruimte en het eiwit niet verder gevouwen kan worden begint de eiwitvouwing helemaal opnieuw. De meest stabiel gevouwen configuratie wordt steeds bewaard. Bekijk [hier](https://github.com/SanneHoeken/The-Triforce-of-Pow-d-er/blob/master/code/algorithms/random.py) de code voor implementatiedetails.
2. **Greedy.**
  Het Greedy algoritme vouwt voor elke iteratie een eiwit door aan elk aminozuur één voor één vouwwaarden toe te wijzen die op dat punt de meest stabiele geldige configuratie oplevert. Als een aminozuur eindigt in een doodlopende ruimte en het eiwit niet verder gevouwen kan worden begint de eiwitvouwing helemaal opnieuw. De meest stabiel gevouwen eindconfiguratie wordt steeds bewaard. Bekijk [hier](https://github.com/SanneHoeken/The-Triforce-of-Pow-d-er/blob/master/code/algorithms/greedy.py) de code voor implementatiedetails.
3. **Hillclimber.**
  Het Hillclimber algoritme neemt een reeds gevouwen eiwit als uitgangspunt en muteert het eiwit een x aantal keer door voor y willekeurige amino's de vouwwaarden te veranderen naar willekeurig gekozen andere waarden. Elke mutatie(serie) die tot een stabielere geldige configuratie leidt wordt bewaard. Bekijk [hier](https://github.com/SanneHoeken/The-Triforce-of-Pow-d-er/blob/master/code/algorithms/hill_climber.py) de code voor implementatiedetails.
4. **Simulated Annealing.**
  Het Simulated Annealing algoritme neemt een reeds gevouwen eiwit als uitgangspunt en muteert het eiwit een x aantal keer door voor y willekeurige amino's de vouwwaarden te veranderen naar willekeurig gekozen andere waarden. Elke mutatie(serie) die tot een stabielere geldige configuratie leidt wordt bewaard. Het Simulated Annealing algoritme accepteert soms ook gemuteerde configuraties die minder stabiel zijn, afhankelijk van de temperatuur die ingesteld is. De temperatuur koelt lineair af met het aantal iteraties lineair. Bekijk [hier](https://github.com/SanneHoeken/The-Triforce-of-Pow-d-er/blob/master/code/algorithms/simulated_annealing.py) de code voor implementatiedetails.
5. **Breadth First Search ++++.**
  Het Breadth First Search algoritme bouwt een boom van de verschillende mogelijkheden om het eiwit te vouwen en doorzoekt deze boom tegelijkertijd in de breedte totdat het de beste vouwmogelijkheid heeft gevonden.
  Een totale opbouw en doorzoeking van de boom zou veel te veel geheugen en tijd vragen, dus dit algoritme is verbeterd met verschillende soorten pruning:
    * Queue restriction: om te voorkomen dat er te veel onnodige mogelijkheden worden bewaard, kan de gebruiker de queue van nodes beperken tot een bepaalde lengte. Dit zal BFS dwingen om slechts een bepaald aantal nodes per rij te behouden en er moet dus verstandig worden gebalanceerd met de andere pruningen, om te voorkomen dat slechts een klein deel van de takken aan de linkerkant van de boom wordt gebruikt.
    * Pruning op relevantie: voor elke mogelijkheid van vouwing wordt een relevantie-score berekend (op basis van de diepte van de node, de huidige score van de aminozuurketen en de structuur van het eiwit op dat punt in de boom). Nodes onder deze score worden gepruned.
    * Pruning op diepte: met een diepteparameter kan de gebruiker een diepte selecteren waaruit nodes op breedte worden gepruned.
    * Pruning op breedte: de gebruiker kan een afstandsfactor geven om meer nodes in de breedte te prunen (om te voorkomen dat er te veel nodes zijn met een vergelijkbare structuur en dezelfde score). Per diepte wordt er slechts één node bewaard iedere diepte * afstandsfactor node. 
  Met gebruik en goede balancering van deze parameters kan het algoritme snel de relevante vouwingen selecteren en een eiwit met een hoge score opbouwen. Bekijk [hier](https://github.com/SanneHoeken/The-Triforce-of-Pow-d-er/blob/master/code/algorithms/bfs_plus.py) code voor implementatiedetails.
6. **Branch & Bound.**
  Het Branch & Bound algoritme is gebaseerd op het werk van Mao Chen & Wen-Qi Huang, 2005. Het Branch & Bound algoritme vouwt een eiwit amino voor amino door te zoeken naar de beste mogelijkheden. In de zoekmethode wordt voor elke H- en C-amino de potentie van de partiële configuratie beoordeeld. Een negatieve evaluatie leidt tot pruning van die configuratie. Bij P-aminos wordt er niet gepruned. 
  De potentie van een configuratie wordt bepaald aan de hand van twee variabelen: 
    * De gemiddelde score van een configuratie met een bepaalde lengte. 
    * De beste score voor een configuratie met een bepaalde lengte. 
  De score van een configuratie wordt vergeleken met deze twee variabelen. Als de score beter is dan de beste score, wordt de partiële configuratie niet gepruned. Een score slechter dan de gemiddelde score wordt gepruned met een gegeven probability p1. Een score beter dan de gemiddelde waarde maar slechter dan de beste score wordt gepruned met een probability p2. Bekijk [hier](https://github.com/SanneHoeken/The-Triforce-of-Pow-d-er/blob/master/code/algorithms/branch_and_bound.py) de code voor implementatiedetails en de bronvermelding.

## Auteurs
- Thomas van Genderen (studentnummer: 11218290)
- Charlotte Lafage (studentnummer: 12977772)
- Sanne Hoeken (studentnummer: 12901202)
