#import "@local/gost732-2017:0.3.0": *


= Анализ требований
#v(1em)

Согласно требованиям технического задания, результатом работы устройства является вычисление сопротивления в диапазоне от 1Ом до 1МОм. Устройство должно быть реализовано на доступной элементной базе с использованием стандартных цифровых компонентов.

== Описание принципа работы описываемого устройства
#v(1em)

В качестве входных данных устройство принимает измеряемый резистор, подключаемый к специальным контактам. Данный резистор используется как элемент R в составе времязадающей RC-цепи, совместно с конденсатором известной емкости C, для получения характерной постоянной времени T = R*C.

Принцип действия разрабатываемого устройства заключается в измерении временного интервала, необходимого для заряда конденсатора C через измеряемое сопротивление R до порогового напряжения. Данный временной интервал напрямую связан с постоянной времени RC-цепи и, следовательно, пропорционален измеряемому сопротивлению R. Постоянная времени определяется как время, за которое напряжение на конденсаторе при его заряде достигает 63.2% от напряжения питания.

Для точного измерения указанного временного интервала используются цифровые счетчики, тактируемые импульсами генератора тактовой частоты. Полученное количество тактовых импульсов за время заряда отображается на цифровом индикаторе как измеренное значение сопротивления, и полученная величина кратна измеряемому сопротивлению.


== Выбор схемотехнического решения
#v(1em)

Исходя из описанного принципа измерения сопротивления, основанного на заряде RC-цепи, выбирается структурная схема устройства, представленная в Приложении Б. 
Процесс измерения инициируется блоком управления, который координирует взаимодействие всех узлов устройства.


По команде начала измерения блок управления формирует сигнал, запускающий заряд конденсатора в измерительной RC-цепи, входящей в состав измерительного узла. Одновременно с началом заряда блок управления разрешает работу блока счетчиков. Блок счетчиков начинает подсчет импульсов, поступающих от генератора тактовых импульсов, который обеспечивает стабильную временную базу для измерений. 


В течение процесса заряда напряжение на конденсаторе непрерывно контролируется блоком сравнения, также являющимся частью измерительного узла. Как только напряжение на конденсаторе достигает заранее установленного порогового значения, блок сравнения формирует сигнал окончания заряда. Этот сигнал немедленно поступает в блок управления. 


Получив сигнал окончания заряда, блок управления формирует сигнал запрета счета, который останавливает работу блока счетчиков. В блоке счетчиков фиксируется итоговое количество тактовых импульсов, подсчитанное строго за время заряда конденсатора до порогового уровня.


Зафиксированное в блоке счетчиков цифровое значение пропорционально измеряемому сопротивлению. Блок управления передает сигнал готовности данных в блок вывода результата. Блок вывода считывает цифровой код из блока счетчиков формирует управляющие сигналы для отображения результата на цифровом индикаторе. 
