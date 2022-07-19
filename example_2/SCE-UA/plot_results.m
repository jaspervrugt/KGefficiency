function plot_results(subset,bestSets,Measurement,Extra,k)

Extra.calPeriod = subset(k,:);
[SimRR] = hymod(bestSets(k,:),Extra);

figure
a = [Extra.calPeriod(1):Extra.calPeriod(2)]';
b = Extra.PET(a);
c = Extra.Precip(a);
subplot(3,1,1)
plot(Extra.numTime(a),b,'-r')
 ylabel('PET')
datetick('x',1)

subplot(3,1,2)
plot(Extra.numTime(a),c,'-c')
ylabel('Precip')
datetick('x',1)

subplot(3,1,3)
plot(Extra.numTime(a),Measurement.MeasData(a),'.',Extra.numTime(a),SimRR)
ylabel('Discharge')
datetick('x',1)