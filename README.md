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
    - **/code/algorithms/help_methods**: bevat benodigde functies voor de algoritmes
  - **/code/classes**: bevat de benodigde classes voor deze case
  - **/code/visualisation**: bevat de code voor de visualisatie
- **/data**: bevat een tekstbestand met eiwit-strings die als input zijn te gebruiken

### Algoritmes

1. Random
  Het Random algoritme vouwt voor elke iteratie een eiwit door willekeurige mogelijke vouwwaarden één voor één toe te wijzen aan elk aminozuur. Als een aminozuur eindigt in een doodlopende ruimte en het eiwit niet verder gevouwen kan worden begint de eiwitvouwing helemaal opnieuw. De meest stabiel gevouwen configuratie wordt steeds bewaard. Bekijk [hier](https://github.com/SanneHoeken/The-Triforce-of-Pow-d-er/blob/master/code/algorithms/random.py) de code voor implementatiedetails.
2. Greedy
  Het Greedy algoritme vouwt voor elke iteratie een eiwit door aan elk aminozuur één voor één vouwwaarden toe te wijzen die op dat punt de meest stabiele geldige configuratie oplevert. Als een aminozuur eindigt in een doodlopende ruimte en het eiwit niet verder gevouwen kan worden begint de eiwitvouwing helemaal opnieuw. De meest stabiel gevouwen eindconfiguratie wordt steeds bewaard. Bekijk [hier](https://github.com/SanneHoeken/The-Triforce-of-Pow-d-er/blob/master/code/algorithms/greedy.py) de code voor implementatiedetails.
3. Hillclimber
  Het Hillclimber algoritme neemt een reeds gevouwen eiwit als uitgangspunt en muteert het eiwit een x aantal keer door voor y willekeurige amino's de vouwwaarden te veranderen naar willekeurig gekozen andere waarden. Elke mutatie(serie) die tot een stabielere geldige configuratie leidt wordt bewaard. Bekijk [hier](https://github.com/SanneHoeken/The-Triforce-of-Pow-d-er/blob/master/code/algorithms/hill_climber.py) de code voor implementatiedetails.
4. Simulated Annealing
  Het Simulated Annealing algoritme neemt een reeds gevouwen eiwit als uitgangspunt en muteert het eiwit een x aantal keer door voor y willekeurige amino's de vouwwaarden te veranderen naar willekeurig gekozen andere waarden. Elke mutatie(serie) die tot een stabielere geldige configuratie leidt wordt bewaard. Het Simulated Annealing algoritme accepteert soms ook gemuteerde configuraties die minder stabiel zijn, afhankelijk van de temperatuur die ingesteld is. De temperatuur koelt met aantal iteraties lineair af. Bekijk [hier](https://github.com/SanneHoeken/The-Triforce-of-Pow-d-er/blob/master/code/algorithms/simulated_annealing.py) de code voor implementatiedetails.
5. Breadth First Search ++++
  Hier informatie over het algoritme. Bekijk [hier](link naar de code) code voor implementatiedetails.
6. Depth First Search met Branch & Bound
  Hier informatie over het algoritme. Bekijk [hier](link naar de code) de code voor implementatiedetails.

## Auteurs
- Thomas van Genderen
- Charlotte Lafage
- Sanne Hoeken