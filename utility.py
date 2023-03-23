from os import path
from enum import Enum

class ColorCode(Enum):
    BLACK = 'black'
    RED = 'salmon'
    GREEN = 'limegreen'
    BLUE = 'cornflowerblue'
    YELLOW = 'yellow'
    PURPLE = 'violet'
    GRAY = 'darkgray'
    
similar_color_dict = dict({ 
    ColorCode.BLACK:('black', 'dimgray'),
    ColorCode.RED:('red','salmon'),
    ColorCode.GREEN:('green','lightgreen'),
    ColorCode.BLUE:('blue', 'cornflowerblue'),
    ColorCode.YELLOW:('gold', 'yellow' ),
    ColorCode.PURPLE:('purple', 'violet' ),
})

def ToOrdinal(n):
    if n//10 !=1 and 1<=n%10<=3  :
        return str(n)+["st","nd","rd"][n-1]
    else:
        return str(n)+"th"

def RelativePathToAbs( p ):
    if path.isabs(p) :
        return p
    return path.abspath(p)
    
def UnpackListOfTuple( list_of_tuple ):
    return [list(tup) for tup in zip(*list_of_tuple)]