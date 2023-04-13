VALUE2BAND  = { 0 : 'water',
                1 : 'trees',
                2 : 'grass',
                3 : 'flooded_vegetation',
                4 : 'crops',
                5 : 'shrub_and_scrub',
                6 : 'built',
                7 : 'bare',
                8 : 'snow_and_ice' }

BAND2VALUE  = { 'water' : 0,
                'trees' : 1,
                'grass' : 2,
                'flooded_vegetation' : 3,
                'crops' : 4,
                'shrub_and_scrub' : 5,
                'built' : 6,
                'bare' : 7,
                'snow_and_ice' : 8 }

                
ELEGIBLE_BANDS = ['grass', 'crops', 'shrub_and_scrub', 'bare']

HECT_CONVERTION_FACTOR = .001

CARBON_RATE = 14.55

MAX_RES_SCALE = 10 # 1 pixel = 10x10m